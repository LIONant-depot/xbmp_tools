namespace xbmp::tools::loader
{
    enum class error : std::uint8_t
    {   SUCCESS
    ,   FAILURE
    };

    error* LoadDSS      ( xcore::bitmap& Bitmap, const char* pFileName ) noexcept;
    error* LoadDSS      ( xcore::bitmap& Bitmap, std::span<const std::byte> Buffer ) noexcept;
    error* LoadSTDImage ( xcore::bitmap& Bitmap, const char* pFileName ) noexcept;
}