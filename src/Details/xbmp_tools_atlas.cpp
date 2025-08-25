namespace xbmp::tools{

//-------------------------------------------------------------------------------------------------------
// int base on size
//-------------------------------------------------------------------------------------------------------
template< std::size_t T_SIZE_BYTES >
using byte_size_uint_t = std::tuple_element_t< T_SIZE_BYTES - 1, std::tuple<std::uint8_t, std::uint16_t, std::uint32_t, std::uint32_t, std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t>>;

template< std::size_t T_SIZE_BYTES >
using byte_size_int_t = std::tuple_element_t< T_SIZE_BYTES - 1, std::tuple<std::int8_t, std::int16_t, std::int32_t, std::int32_t, std::int64_t, std::int64_t, std::int64_t, std::int64_t>>;

//-------------------------------------------------------------------------------------------------------
// given a err returns the equivalent size as either a sing or unsigned err
//-------------------------------------------------------------------------------------------------------
template< typename T >
using to_int_t = byte_size_int_t<sizeof(T)>;

template< typename T >
using to_uint_t = byte_size_uint_t<sizeof(T)>;

//------------------------------------------------------------------------------
// Description:
//      Takes an Address or integer and aligns it up base on the given alignment.
//      The result in the next number greater than or equal to "n" 
//      which is a multiple of "a".  For example, x_Align(57,16) is 64.
// Arguments:
//     Addr       - This is the address/number/offset to align
//     AlignTo    - This is a power of 2 number that the user_types wants it to be align to
//------------------------------------------------------------------------------
template< typename T > constexpr 
T Align(T Address, const int AlignTo) noexcept
{
    static_assert(std::is_integral<T>::value, "This function only works with integer values");
    using unsigned_t = to_uint_t<T>;
    return static_cast<T>((unsigned_t(Address) + (static_cast<unsigned_t>(AlignTo) - 1)) & static_cast<unsigned_t>(-AlignTo));
}


//------------------------------------------------------------------------------
template< typename T > constexpr
T NextPowOfTwo(const T x, const int s) noexcept
{
    static_assert(std::is_integral<T>::value, "");
    return static_cast<T>((s == 0) ? 1 + x : NextPowOfTwo<T>(x | (x >> s), s >> 1));
}

//------------------------------------------------------------------------------
// Description:
//      Rounds a number to the next power of two
// Example:
//      RoundToNextPowOfTwo(3) == 4   // The next power from #3 is #4
//------------------------------------------------------------------------------
template< typename T> constexpr
T RoundToNextPowOfTwo(const T x) noexcept
{
    static_assert(std::is_integral<T>::value, "");
    return (x == 0) ? 0 : NextPowOfTwo<T>(x - 1, static_cast<int>(4 * sizeof(T)));
}


//-------------------------------------------------------------------------------------------------
static
int Area( const atlas::rect_xywhf** pA, const atlas::rect_xywhf** pB )
{
    int b1 = (*pA)->getArea();
    int a1 = (*pB)->getArea();
    if( a1 < b1 ) return -1;
    return a1>b1;
}

//-------------------------------------------------------------------------------------------------
static
int Perimeter( const atlas::rect_xywhf** pA, const atlas::rect_xywhf** pB )
{
    int b1 = (*pA)->getPerimeter();
    int a1 = (*pB)->getPerimeter();
    if( a1 < b1 ) return -1;
    return a1>b1;
}

//-------------------------------------------------------------------------------------------------
static
int MaxSide( const atlas::rect_xywhf** pA, const atlas::rect_xywhf** pB )
{
    int b1 = std::max( (*pA)->m_W, (*pA)->m_H );
    int a1 = std::max( (*pB)->m_W, (*pB)->m_H );
    if( a1 < b1 ) return -1;
    return a1>b1;
}

//-------------------------------------------------------------------------------------------------
static
int MaxWidth( const atlas::rect_xywhf** pA, const atlas::rect_xywhf** pB )
{
    int b1 = (*pA)->m_W;
    int a1 = (*pB)->m_W;
    if( a1 < b1 ) return -1;
    return a1>b1;
}

//-------------------------------------------------------------------------------------------------
static
int MaxHeight( const atlas::rect_xywhf** pA, const atlas::rect_xywhf** pB )
{
    int b1 = (*pA)->m_H;
    int a1 = (*pB)->m_H;
    if( a1 < b1 ) return -1;
    return a1>b1;
}

// just add another comparing function name to cmpf to perform another packing attempt
// more functions == slower but probably more efficient cases covered and hence less area wasted
constexpr static auto comparison_fn_v = std::array
{
    Area,
    Perimeter,
    MaxSide,
    MaxWidth,
    MaxHeight
};

// if you find the algorithm running too slow you may double this factor to increase speed but also decrease efficiency
// 1 == most efficient, slowest
// efficiency may be still satisfying at 64 or even 256 with nice speedup

constexpr static int s_discard_step_v = 16;

// For every sorting function, algorithm will perform packing attempts beginning with a bin with width
// and height equal to max_side, and decreasing its dimensions if it finds out that rectangles did
// actually fit, increasing otherwise. Although, it's doing that in sort of binary search manner, so
// for every comparing function it will perform at most log2(max_side) packing attempts looking for
// the smallest possible bin size. discard_step means that if the algorithm will break of the searching
// loop if the rectangles fit but "it may be possible to fit them in a bin smaller by discard_step"
//
//  may be pretty slow in debug mode anyway (std::vector and stuff like that in debug mode is always slow)
//
// the algorithm was based on http://www.blackpawn.com/texts/lightmaps/default.html
// the algorithm reuses the node tree so it doesn't reallocate them between searching attempts
//
//
// please let me know about bugs at  <----
// unknownunreleased@gmail.com           |
// |
// |
// ps. I'm 16 so take this --------------- more than seriously, though I made some realtime tests
// with packing several hundreds of rectangles every frame, no crashes, no memory leaks, good results
// Thank you.


struct node
{
    struct pnode
    {
        node*   m_pNode;
        bool    m_Fill;
        
        pnode() : m_Fill(false), m_pNode( nullptr ) {}
        
        void setup( int L, int T, int R, int B )
        {
            if( m_pNode )
            {
                m_pNode->m_RC  = atlas::rect_ltrb( L, T, R, B );
                m_pNode->m_pID = nullptr;
            }
            else
            {
                m_pNode = new node();
                m_pNode->setup( atlas::rect_ltrb( L, T, R, B ) );
            }
            
            m_Fill = true;
        }
    };
    
    pnode                 m_C[2];
    atlas::rect_ltrb      m_RC;
    atlas::rect_xywhf*    m_pID;
    
    ~node()
    {
        if( m_C[0].m_pNode ) delete m_C[0].m_pNode;
        if( m_C[1].m_pNode ) delete m_C[1].m_pNode;
    }
    
    node( void )
    {
        setup( atlas::rect_ltrb(0,0,0,0) );
    }
    
    void setup( const atlas::rect_ltrb& RC )
    {
        m_pID = nullptr;
        m_RC  = RC;
    }
    
    void Reset( const atlas::rect_wh& R )
    {
        m_pID = nullptr;
        m_RC  = atlas::rect_ltrb( 0, 0, R.m_W, R.m_H );
        DelCheck();
    }
    
    void DelCheck( void )
    {
        if( m_C[0].m_pNode ) { m_C[0].m_Fill = false; m_C[0].m_pNode->DelCheck(); }
        if( m_C[1].m_pNode ) { m_C[1].m_Fill = false; m_C[1].m_pNode->DelCheck(); }
    }
    
    node* Insert( atlas::rect_xywhf& ImageRect )
    {
        if( m_C[0].m_pNode && m_C[0].m_Fill )
        {
            node* pNewNode = m_C[0].m_pNode->Insert( ImageRect );
            if( pNewNode ) return pNewNode;
            
            pNewNode = m_C[1].m_pNode->Insert( ImageRect );
            return pNewNode;
        }
        
        if( m_pID ) return nullptr;
        
        int f = ImageRect.Fits( atlas::rect_xywh(m_RC) );
        
        switch( f )
        {
            case 0: return 0;
            case 1: ImageRect.m_bFlipped = false; break;
            case 2: ImageRect.m_bFlipped = true;  break;
            case 3: m_pID = &ImageRect; ImageRect.m_bFlipped = false; return this;
            case 4: m_pID = &ImageRect; ImageRect.m_bFlipped = true;  return this;
        }
        
        int iW = (ImageRect.m_bFlipped ? ImageRect.m_H : ImageRect.m_W);
        int iH = (ImageRect.m_bFlipped ? ImageRect.m_W : ImageRect.m_H);
        
        if( ( m_RC.getWidth() - iW) >  ( m_RC.getHeight() - iH) )
        {
            m_C[0].setup( m_RC.m_Left,      m_RC.m_Top,     m_RC.m_Left + iW,   m_RC.m_Bottom );
            m_C[1].setup( m_RC.m_Left + iW, m_RC.m_Top,     m_RC.m_Right,       m_RC.m_Bottom );
        }
        else
        {
            m_C[0].setup( m_RC.m_Left,      m_RC.m_Top,      m_RC.m_Right,      m_RC.m_Top + iH   );
            m_C[1].setup( m_RC.m_Left,      m_RC.m_Top + iH, m_RC.m_Right,      m_RC.m_Bottom     );
        }
        
        return m_C[0].m_pNode->Insert( ImageRect );
    }
};

//-------------------------------------------------------------------------------------------------
static
atlas::rect_wh TheRect2D( atlas::rect_xywhf* const*             pV
                        , atlas::pack_mode                      Mode
                        , int                                   n
                        , int                                   MaxS
                        , std::vector<atlas::rect_xywhf*>&      Succ
                        , std::vector<atlas::rect_xywhf*>&      UnSucc )
{
    node                    Root;
    const int               nFuncs = 5;
    atlas::rect_xywhf**     Order[ nFuncs ];
    const int               Multiple = (Mode == atlas::pack_mode::MULTIPLE_OF_4)?4:8;
    
    for( int f = 0; f < nFuncs; ++f )
    {
        Order[f] = new atlas::rect_xywhf*[n]{};
        std::memcpy( Order[f], pV, sizeof(atlas::rect_xywhf*) * n);
        std::qsort( (void*)Order[f], n, sizeof(atlas::rect_xywhf*), (int(*)(const void*, const void*))comparison_fn_v[f] );
    }
    
    atlas::rect_wh    MinBin      = atlas::rect_wh( MaxS, MaxS );
    int                     MinFunc     = -1;
    int                     BestFunc    = 0;
    int                     BestArea    = 0;
    int                     TheArea     = 0;
    bool                    bFail       = false;
    int                     Step, Fit;
    
    for( int f = 0; f < nFuncs; ++f )
    {
        pV      = Order[f];
        Step    = MinBin.m_W / 2;
        Root.Reset( MinBin );
        
        while(1)
        {
            if( Root.m_RC.getWidth() > MinBin.m_W )
            {
                if( MinFunc > -1 )
                    break;
                TheArea = 0;
                
                Root.Reset( MinBin );
                
                for( int i = 0; i < n; ++i )
                {
                    if( Root.Insert( *pV[i]) )
                        TheArea += pV[i]->getArea();
                }
                
                bFail = true;
                break;
            }
            
            Fit = -1;
            
            for( int i = 0; i < n; ++i )
            {
                if( !Root.Insert( *pV[i]) )
                {
                    Fit = 1;
                    break;
                }
            }
            
            if( Fit == -1 && Step <= s_discard_step_v )
                break;
            
            int OldW = Root.m_RC.getWidth();
            int OldH = Root.m_RC.getHeight();
            int NewW = OldW  + Fit*Step;
            int NewH = OldH  + Fit*Step;
            
            if( Mode == atlas::pack_mode::MULTIPLE_OF_4 ||
                Mode == atlas::pack_mode::MULTIPLE_OF_8 )
            {
                NewW = Align( NewW, Multiple );
                NewH = Align( NewH, Multiple );
                
                assert( Align( OldW, Multiple) == OldW );
                assert( Align( OldH, Multiple) == OldH );
            }
            else if( Mode == atlas::pack_mode::POWER_OF_TWO ||
                     Mode == atlas::pack_mode::POWER_OF_TWO_SQUARE )
            {
                NewW = RoundToNextPowOfTwo( NewW );
                NewH = RoundToNextPowOfTwo( NewH );
                
                assert( RoundToNextPowOfTwo( OldW ) == OldW );
                assert( RoundToNextPowOfTwo( OldH ) == OldH );
            }

            
            //
            // Try to choose the search in the right dimension
            //
            if( Mode != atlas::pack_mode::ANY_SIZE )
            {
                if( Fit == 1 )
                {
                    if( OldW < OldH ) NewH = OldH;
                    else              NewW = OldW;
                }
                else
                {
                    if( OldW < OldH ) NewW = OldW;
                    else              NewH = OldH;
                }
            }
            
            //
            // For square sizes we need to pick the best dimension for all
            //
            if( Mode == atlas::pack_mode::POWER_OF_TWO_SQUARE )
            {
                if( Fit == 1 )
                {
                    NewH = NewW = std::max( NewW, NewH );
                }
                else
                {
                    NewH = NewW = std::min( NewW, NewH );
                }
            }

            //
            // Ok now we can change the size
            //
            Root.Reset( atlas::rect_wh( NewW, NewH ) );

            //
            // Keep searching
            //
            Step /= 2;
            if( !Step )
                Step = 1;
        }
        
        if( !bFail && ( MinBin.getArea() >= Root.m_RC.getArea()) )
        {
            MinBin  = atlas::rect_wh( Root.m_RC );
            MinFunc = f;
        }
        else if( bFail && ( TheArea > BestArea) )
        {
            BestArea = TheArea;
            BestFunc = f;
        }
        
        bFail = false;
    }
    
    pV = Order[ MinFunc == -1 ? BestFunc : MinFunc ];
    
    int     ClipX = 0, ClipY = 0;
    node*   pRet;
    
    Root.Reset( MinBin );
    
    for( int i = 0; i < n; ++i )
    {
        pRet = Root.Insert( *pV[i]);
        if( pRet )
        {
            pV[i]->m_X = pRet->m_RC.m_Left;
            pV[i]->m_Y = pRet->m_RC.m_Top;
            
            if( pV[i]->m_bFlipped )
            {
                pV[i]->m_bFlipped = false;
                pV[i]->Flip();
            }
            
            ClipX = std::max( ClipX, pRet->m_RC.m_Right );
            ClipY = std::max( ClipY, pRet->m_RC.m_Bottom );
            
            Succ.push_back(pV[i]);
        }
        else
        {
            UnSucc.push_back(pV[i]);
            
            pV[i]->m_bFlipped = false;
        }
    }
    
    for( int f = 0; f < nFuncs; ++f )
    {
        delete Order[f];
    }
    
    // Make sure that we do return a power of two texture
    if( Mode == atlas::pack_mode::POWER_OF_TWO ||
        Mode == atlas::pack_mode::POWER_OF_TWO_SQUARE )
    {
        return atlas::rect_wh( RoundToNextPowOfTwo(ClipX), RoundToNextPowOfTwo(ClipY) );
    }
    else if( Mode == atlas::pack_mode::MULTIPLE_OF_4 ||
             Mode == atlas::pack_mode::MULTIPLE_OF_8 )
    {
        const int Multiple = (Mode == atlas::pack_mode::MULTIPLE_OF_4)?4:8;
        return atlas::rect_wh( Align(ClipX, Multiple), Align(ClipY,Multiple) );
    }
    
    return atlas::rect_wh( (ClipX), (ClipY) );
}

//-------------------------------------------------------------------------------------------------

bool atlas::Pack( rect_xywhf* const* pV, int nRects, int MaxSize, std::vector<bin>& Bins, pack_mode Mode )
{
    rect_wh TheRect( MaxSize, MaxSize );
    
    // make sure that all the rectangles fit in the given MaxSize
    for( int i = 0; i < nRects; ++i )
    {
        if( !pV[i]->Fits(TheRect) )
            return false;
    }
    
    // Create a double buffer array of pointers
    std::vector<rect_xywhf*>     Vec[2];
    std::vector<rect_xywhf*>*    pVec[2] = { &Vec[0], &Vec[1] };
    
    // Set one pointer per rectangle past by the user
    Vec[0].resize( nRects );
    
    // Initialize the firt array to zero
    std::memcpy( &Vec[0][0], pV, sizeof(rect_xywhf*) * nRects );
    
    while(1)
    {
        bin& NewBin = Bins.emplace_back();
        
        NewBin.m_Size = TheRect2D( &((*pVec[0])[0]), Mode, static_cast<int>(pVec[0]->size()), MaxSize, NewBin.m_Rects, *pVec[1] );
        
        pVec[0]->clear();
        
        if( pVec[1]->size() == 0 ) break;
        
        std::swap( pVec[0], pVec[1] );
    }
    
    return true;
}

//-------------------------------------------------------------------------------------------------

atlas::rect_wh::rect_wh(const rect_ltrb& rr) : m_W(rr.getWidth()), m_H(rr.getHeight()) {}
atlas::rect_wh::rect_wh(const rect_xywh& rr) : m_W(rr.getWidth()), m_H(rr.getHeight()) {}
atlas::rect_wh::rect_wh(int w, int h) : m_W(w), m_H(h) {}

//-------------------------------------------------------------------------------------------------

int atlas::rect_wh::Fits( const rect_wh& r ) const
{
    if( m_W == r.m_W && m_H == r.m_H ) return 3;
    if( m_H == r.m_W && m_W == r.m_H ) return 4;
    if( m_W <= r.m_W && m_H <= r.m_H ) return 1;
    if( m_H <= r.m_W && m_W <= r.m_H ) return 2;
    return 0;
}

//-------------------------------------------------------------------------------------------------

atlas::rect_ltrb::rect_ltrb( void ) : m_Left(0), m_Top(0), m_Right(0), m_Bottom(0) {}
atlas::rect_ltrb::rect_ltrb( int l, int t, int r, int b) : m_Left(l), m_Top(t), m_Right(r), m_Bottom(b) {}

//-------------------------------------------------------------------------------------------------

int atlas::rect_ltrb::getWidth( void ) const
{
    return m_Right - m_Left;
}

//-------------------------------------------------------------------------------------------------

int atlas::rect_ltrb::getHeight( void ) const
{
    return m_Bottom - m_Top;
}

//-------------------------------------------------------------------------------------------------

int atlas::rect_ltrb::getArea( void ) const
{
    return getWidth() * getHeight();
}

//-------------------------------------------------------------------------------------------------

int atlas::rect_ltrb::getPerimeter( void ) const
{
    return 2*getWidth() + 2*getHeight();
}

//-------------------------------------------------------------------------------------------------

void atlas::rect_ltrb::setWidth( int Width )
{
    m_Right = m_Left + Width;
}

//-------------------------------------------------------------------------------------------------

void atlas::rect_ltrb::setHeight( int Height )
{
    m_Bottom = m_Top + Height;
}

//-------------------------------------------------------------------------------------------------

atlas::rect_xywh::rect_xywh( void ) : m_X(0), m_Y(0) {}
atlas::rect_xywh::rect_xywh( const rect_ltrb& rc ) :
    m_X( rc.m_Left ), m_Y( rc.m_Top )
{
    m_H = rc.m_Bottom - rc.m_Top;
    m_W = rc.m_Right  - rc.m_Left;
}

//-------------------------------------------------------------------------------------------------

atlas::rect_xywh::rect_xywh( int x, int y, int w, int h ) : m_X(x), m_Y(y), rect_wh( w, h ) {}

//-------------------------------------------------------------------------------------------------

atlas::rect_xywh::operator atlas::rect_ltrb( void )
{
    rect_ltrb rr( m_X, m_Y, 0, 0 );
    rr.setWidth ( m_W );
    rr.setHeight( m_H );
    return rr;
}

//-------------------------------------------------------------------------------------------------

int atlas::rect_xywh::getRight( void ) const
{
    return m_X + m_W;
};

//-------------------------------------------------------------------------------------------------

int atlas::rect_xywh::getBottom( void ) const
{
    return m_Y + m_H;
}

//-------------------------------------------------------------------------------------------------

void atlas::rect_xywh::setWidth( int Right )
{
    m_W = Right - m_X;
}

//-------------------------------------------------------------------------------------------------

void atlas::rect_xywh::setHeight( int Bottom )
{
    m_H = Bottom - m_Y;
}

//-------------------------------------------------------------------------------------------------

int atlas::rect_wh::getArea( void ) const
{
    return m_W * m_H;
}

//-------------------------------------------------------------------------------------------------

int atlas::rect_wh::getPerimeter( void ) const
{
    return 2 * m_W + 2 * m_H;
}

//-------------------------------------------------------------------------------------------------

atlas::rect_xywhf::rect_xywhf( const rect_ltrb& rr ) : rect_xywh(rr), m_bFlipped( false ) {}
atlas::rect_xywhf::rect_xywhf( int x, int y, int width, int height ) : rect_xywh( x, y, width, height ), m_bFlipped( false ) {}
atlas::rect_xywhf::rect_xywhf() : m_bFlipped( false ) {}

//-------------------------------------------------------------------------------------------------

void atlas::rect_xywhf::Flip( void )
{
    m_bFlipped = !m_bFlipped;
    std::swap( m_W, m_H );
}

} // end of namespace xbmp::tools 