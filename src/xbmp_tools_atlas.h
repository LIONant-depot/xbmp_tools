namespace xbmp::tools
{
    // rect_xywhf - structure representing your rectangle object
    // bin - structure representing your bin object
    // bool rect2D(rect_xywhf* const * v, int n, int max_side, std::vector<bin>& bins) - actual packing function
    //
    // v        - pointer to array of pointers to your rectangle (" ** " design is for the
    //            sake of sorting speed and const means that the pointers will point to the same rectangles)
    // n        - pointers' count
    // max_side - maximum bins' side - algorithm works with square bins (in the end it may trim them
    //            to rectangular form). for the algorithm to finish faster, pass a reasonable value
    //            (unreasonable would be passing 1 000 000 000 for packing 4 50x50 rectangles).
    // bins     - vector to which the function will push_back() created bins, each of them containing vector
    //            to pointers of rectangles from "v" belonging to that particular bin.
    //            Every bin also keeps information about its width and height of course, none of the
    //            dimensions is bigger than max_side.
    //
    // returns true on success, false if one of the rectangles' dimension was bigger than max_side
    //
    // You want to pass your vector of rectangles representing your textures/glyph objects with
    // GL_MAX_TEXTURE_SIZE as max_side, then for each bin iterate through its rectangles, typecast each
    // one to your own structure and then memcpy its pixel contents (rotated by 90 degrees if "flipped"
    // rect_xywhf's member is true) to the array representing your texture atlas to the place specified
    // by the rectangle, then in the end upload it with glTexImage2D.
    //
    // Algorithm doesn't create any new rectangles.
    // You just pass an array of pointers and rectangles' x/y/w/h are modified in place, with just vector
    // of pointers for every new bin to let you know which ones belong to the particular bin.
    // Modifying w/h means that the dimensions can be actually swapped for the sake of fitting, the flag
    // "flipped" will be set to true if such a thing occurs.
    //
    // For description how to tune the algorithm and how it actually works see the .cpp file.
    //
    class atlas
    {
    public:

        struct rect_ltrb;
        struct rect_xywh;

        struct rect_wh
        {
                        // 0 - no, 1 - yes, 2 - flipped, 3 - perfectly, 4 perfectly flipped
                        rect_wh         ( const rect_ltrb& );
                        rect_wh         ( const rect_xywh& );
                        rect_wh         ( int W = 0, int H = 0 );
            int         Fits            ( const rect_wh& bigger ) const;
            int         getArea         ( void ) const;
            int         getPerimeter    ( void ) const;
            
            int         m_W;
            int         m_H;
            int         m_Area;
            int         m_Perimeter;
        };

        // rectangle implementing left/top/right/bottom behaviour

        struct rect_ltrb
        {
                        rect_ltrb           ( void );
                        rect_ltrb           ( int left, int top, int right, int bottom );
            
            int         m_Left;
            int         m_Top;
            int         m_Right;
            int         m_Bottom;
            
            int         getWidth            ( void ) const;
            int         getHeight           ( void ) const;
            int         getArea             ( void ) const;
            int         getPerimeter        ( void ) const;
            void        setWidth            ( int Width );
            void        setHeight           ( int Height );
        };

        struct rect_xywh : public rect_wh
        {
                        rect_xywh           ( void );
                        rect_xywh           ( const rect_ltrb& );
                        rect_xywh           ( int x, int y, int width, int height );
            int         getWidth            ( void ) const { return m_W; }
            int         getHeight           ( void ) const { return m_H; }
            int         getRight            ( void ) const;
            int         getBottom           ( void ) const;
            void        setWidth            ( int Right );
            void        setHeight           ( int Bottom );    
                        operator rect_ltrb  ( void );
            
            int         m_X;
            int         m_Y;
        };

        struct rect_xywhf : public rect_xywh
        {
                        rect_xywhf          ( const rect_ltrb& );
                        rect_xywhf          ( int X, int Y, int Width, int Height );
                        rect_xywhf          ( void );
            void        Flip                ( void );
            
            bool        m_bFlipped;
            rect_xywhf* m_pNextFree;
        };

        struct bin
        {
            rect_wh                    m_Size;
            xcore::vector<rect_xywhf*> m_Rects;
        };
        
        enum class pack_mode : std::uint8_t
        { ANY_SIZE
        , MULTIPLE_OF_4
        , MULTIPLE_OF_8
        , POWER_OF_TWO
        , POWER_OF_TWO_SQUARE
        , ENUM_COUNT
        };

        bool  Pack( rect_xywhf* const* pV, int nRects, int MaxSide, xcore::vector<bin>& bins, pack_mode Mode = pack_mode::ANY_SIZE);
    };
}