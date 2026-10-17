#pragma once
//
// File access by name for xbmp_tools.
// File names come in as wide strings (UTF-16 on Windows). On Windows a narrow name given to fopen or std::ifstream is
// read in the ANSI code page, so any character outside it (Cyrillic, CJK, ...) is lost and the file is not found.
// The libraries used here (stb_image, stb_image_write, tinyddsloader) only take narrow names, so the file is opened
// here (_wfopen_s on Windows, UTF-8 elsewhere) and the libraries get the FILE* or the bytes.
//
#include <cstdio>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "dependencies/xstrtool/source/xstrtool.h"

namespace xbmp::tools::details
{
    struct file_close
    {
        void operator()(std::FILE* pFile) const noexcept { std::fclose(pFile); }
    };

    using unique_file = std::unique_ptr<std::FILE, file_close>;

    //------------------------------------------------------------------------------
    // Opens FileName with an fopen mode (L"rb", L"wb", ...); empty if it can not be opened

    inline unique_file OpenFile(std::wstring_view FileName, const wchar_t* pMode) noexcept
    {
    #if defined(_WIN32)
        std::FILE* pFile = nullptr;
        if (_wfopen_s(&pFile, std::wstring(FileName).c_str(), pMode)) return {};
        return unique_file{ pFile };
    #else
        return unique_file{ std::fopen(xstrtool::UTF16ToUTF8(FileName).c_str(), xstrtool::UTF16ToUTF8(pMode).c_str()) };
    #endif
    }

    //------------------------------------------------------------------------------
    // Reads the rest of an open file; false if reading fails

    inline bool ReadFile(std::vector<std::uint8_t>& Data, std::FILE& File) noexcept
    {
        Data.clear();

        std::uint8_t Buffer[64 * 1024];
        while (const auto n = std::fread(Buffer, 1, sizeof(Buffer), &File))
            Data.insert(Data.end(), Buffer, Buffer + n);

        return 0 == std::ferror(&File);
    }

    //------------------------------------------------------------------------------
    // Destination for the stbi_write_*_to_func functions.
    // Like stbi_write_* the file is only created once there is something to write.

    struct file_writer
    {
        std::wstring_view   m_FileName;
        unique_file         m_File      = {};
        bool                m_bFailed   = false;

        static void Write(void* pContext, void* pData, int Size) noexcept
        {
            auto& This = *static_cast<file_writer*>(pContext);
            if (This.m_bFailed) return;

            if (!This.m_File && !(This.m_File = OpenFile(This.m_FileName, L"wb")))
            {
                This.m_bFailed = true;
                return;
            }

            if (std::fwrite(pData, 1, static_cast<std::size_t>(Size), This.m_File.get()) != static_cast<std::size_t>(Size))
                This.m_bFailed = true;
        }

        // true if the file was created, everything was written and it closed fine
        bool Close() noexcept
        {
            bool bOk = m_File && !m_bFailed;
            if (m_File && std::fclose(m_File.release())) bOk = false;
            return bOk;
        }
    };
}
