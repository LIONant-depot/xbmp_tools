namespace xbmp::tools::writers
{
    enum class error : std::uint8_t
    { SUCCESS
    , FAILURE
    };

    error* SaveSTDImage(const char* pFileName, const xcore::bitmap& Bitmap) noexcept;
}