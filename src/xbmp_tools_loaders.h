namespace xbmp::tools::loader
{
    enum class error : std::uint8_t
    {   SUCCESS
    ,   FAILURE
    };

    error* LoadDSS      ( xcore::bitmap& Bitmap, const char* pFileName ) noexcept;
    error* LoadSTDImage ( xcore::bitmap& Bitmap, const char* pFileName ) noexcept;
}