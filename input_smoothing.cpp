#include "input_smoothing.h"
#include <sstream>

namespace {

inline double sqr( double x ) { return x * x; }

inline int sign( double x ) { return x > 0.0 ? 1 : -1; }

inline double three_pt_deriv( const vector<Point2> &points, size_t at ) {
    auto p0 = points[at - 1], p1 = points[at + 1];
    return ( p1.y - p0.y ) / ( p1.x - p0.x );
}

inline double curvature_from_derivs( double deriv1, double deriv2 ) {
    return deriv2 / std::pow( 1 + sqr( deriv1 ), 3.0 / 2.0 );
}

inline double five_pt_deriv( const vector<Point2> &points, size_t at ) {
    auto p_m2 = points[at - 2], p_m1 = points[at - 1], p_p1 = points[at + 1],
         p_p2 = points[at + 2], p = points[at];
    const auto r =
        ( 1.0 / ( p_p1.x - p_m1.x ) ) * ( ( p_p2.y - p.y ) / ( p_p2.x - p.x ) -
                                          ( p.y - p_m2.y ) / ( p.x - p_m2.x ) );
    return r;
}

inline double
five_pt_deriv( Point2 p_m2, Point2 p_m1, Point2 p, Point2 p_p1, Point2 p_p2 ) {
    const auto r =
        ( 1.0 / ( p_p1.x - p_m1.x ) ) * ( ( p_p2.y - p.y ) / ( p_p2.x - p.x ) -
                                          ( p.y - p_m2.y ) / ( p.x - p_m2.x ) );
    return r;
}

inline double three_pt_deriv( Point2 p0, Point2 p, Point2 p1 ) {
    // todo: check why p is not used
    return ( p1.y - p0.y ) / ( p1.x - p0.x );
}

// http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
// See also: http://www.qc.edu.hk/math/Advanced%20Level/Point_to_line.htm
// And this: https://github.com/skyrpex/psimpl/blob/master/psimpl.h`
inline double point_line_distance( Point2 p, Point2 l0, Point2 l1 ) {
    vec2d v{l0, l1};
    vec2d r{p, l0};
    const double d = ( v * r ) / len( v );
    return d;
}

// Please see
// https://ocw.mit.edu/ans7870/18/18.013a/textbook/HTML/chapter15/section04.html
double
calculate_i_curvature( Point2 p0, Point2 p1, Point2 p2, Point2 p3, Point2 p4 ) {
    double dx = ( p3.x - p1.x ) / 2.0;
    double dy = ( p3.y - p1.y ) / 2.0;
    double d2x = ( ( 2 * p2.x ) - p0.x - p4.x ) / 4.0;
    double d2y = ( ( 2 * p2.y ) - p0.y - p4.y ) / 4.0;
    vec2d v{dx, dy}, a{d2x, d2y};
    double v2 = ( v * v );
    double k = len( ( a * v2 - v * ( a * v ) ) / ( v2 * v2 ) );
    return k * sign( dx / dy );
}

const char *to_chararray( TracePointClass klass ) {
    switch ( klass ) {
    case TracePointClass::Inflection:
        return "INFLECTION";
    case TracePointClass::SharpEdge:
        return "SHARP EDGE";
    default:
        return "NORMAL";
    }
}

} // anonymous namespace

InputSmoothing::InputSmoothing( double tollerance ) {}

// ------------------------------------------------------------------------------------

// x - nothing available
// p - first derivative available
// q - second derivative available
// c - curvature available
// k - point class available
// L - last processed point

//           L         N
// . . . . . . . . . . .
//   p p p p p p p p p
//     q q q q q q q
//       c c c c c
//   k k k k k k k k k
void InputSmoothing::on_trace_diff_available( Trace trace ) {
    if ( trace.samples.empty() ) {
        return;
    }
    append_trace( trace );

    size_t points_to_process = m_offset + trace.samples.size();
    for ( size_t index = m_offset; index < points_to_process; ++index ) {
        if ( can_be_calculated_now( index ) ) {
            process_next_point( index, calculate_class( index ) );
            m_offset += 1;
        } else if ( !can_be_calculated_later( index ) ) {
            // cannot be caluclated later, because there are not points before.
            // this might be possible to remove in case we start interpolation
            // points before (only after test cases are written).
            process_next_point( index, TracePointClass::Normal );
            m_offset += 1;
        } else {
            // .. waiting for more input or end trace event.
            // there are enough points before, but there are no points
            // after so that we postpone calculating (we cannot classify point
            // as normal,
            // because it may be lie.)
            break; // todo: is it needed here?
        }
    }
}

void InputSmoothing::on_trace_end() {
    printf( "processing end! offset: %d, total_size: %d\n", (int)m_offset,
            (int)m_points_buffer.size() - 1 );
    for ( size_t index = m_offset; index < m_points_buffer.size() - 1;
          ++index ) {
        process_next_point( index, TracePointClass::Normal );
    }
    process_next_point( m_points_buffer.size() - 1,
                        TracePointClass::SharpEdge );
}

void InputSmoothing::process_next_point( size_t index,
                                         TracePointClass _class ) {
    printf( "%d: %s  '%s'\n", (int)index, to_chararray( _class ),
            _class != TracePointClass::Normal
                ? fetch_curvature_around( index, 1, 1 ).c_str()
                : "" );
    if ( _class != TracePointClass::Normal ) {
        if ( index - m_segment_begin_index == 0 ) {
            // should be first point of segment
            m_segment_begin_index = index;
            m_segment_tangent0 = fetch_point( index );
            printf( "emit: %s\n", to_string( fetch_point( index ) ).c_str() );
        } else {
            // run DP
            printf( "dbg: encountered critical point\n" );
            printf( "emit: %s\n", to_string( fetch_point( index ) ).c_str() );
            // TODO: Should this one be kept?
            m_segment_begin_index = index;
            m_segment_tangent0 = fetch_point( index );
        }
    } else {
        // Normal, line but current point might be outside of corridor
        if ( index - m_segment_begin_index == 1 ) {
            // just second point of current segment
            printf( "dbg: second point of segment: %d\n", (int)index );
            m_segment_tangent1 = fetch_point( index );
            printf( "dbg: got tangent line!: %s..%s\n",
                    to_string( m_segment_tangent0 ).c_str(),
                    to_string( m_segment_tangent1 ).c_str() );
        } else {
        }

        auto p = fetch_point( index );
        double d =
            point_line_distance( p, m_segment_tangent0, m_segment_tangent1 );
        double d2 = d * d;
        printf( "dbg: distance from point %s to line %s..%s is %4.2f\n",
                to_string( p ).c_str(), to_string( m_segment_tangent0 ).c_str(),
                to_string( m_segment_tangent1 ).c_str(), d );
        if ( d2 > sqr( 3. ) ) {
            printf( "dbg: outside!!!\n" );
            printf( "emit: %s\n", to_string( p ).c_str() );
            // corridor has broken
            m_segment_begin_index = index;
            m_segment_tangent0 = p;
        } else {
            printf( "dbg: inside\n" );
        }
    }
}

string InputSmoothing::fetch_curvature_around( size_t around, size_t back_n,
                                               size_t forth_n ) {
    stringstream ss;
    size_t begin = std::max( (int)around - (int)back_n, 0 );
    size_t end = std::min( around + forth_n, m_points_buffer.size() - 1 ) + 1;
    ss << "Curv(" << begin << ".." << ( end - 1 ) << ")=[";
    ss << setprecision( 4 );
    for ( size_t i = begin; i != end; ++i ) {
        ss << m_curvature[i] << ( i == end - 1 ? "]" : ", " );
    }
    return ss.str();
}

Point2 InputSmoothing::fetch_point( size_t at ) const {
    return m_points_buffer.at( at ); // todo: use [] once it is debugged
}

size_t InputSmoothing::points_before_count( size_t index ) const {
    // this may change once we introduce interpolating prior points
    // (corresponding fetch_point() is also changed in this case)
    return index;
}

size_t InputSmoothing::points_after_count( size_t index ) const {
    // this might change when we decide to interpolate points after trace
    // end
    // (but this is unlikely)
    return m_points_buffer.size() - index - 1;
}

bool InputSmoothing::can_be_calculated_later( size_t index ) const {
    assert( !can_be_calculated_now( index ) );
    return points_before_count( index ) >= 2;
}

bool InputSmoothing::can_be_calculated_now( size_t index ) const {
    return ( points_before_count( index ) >= 2 &&
             points_after_count( index ) >= 2 ) ||
           index == 0;
}

// ------------------------------------------------------------------------
// From the diagram below it should be clear, that we can compute
// class only from point 4, because we need smoothed curvature values,
// which need one prior i-curvature value and each curvature needs
// first and second derivative and second derivative requires two points
// before and two after.
// Further improvement could be interpolating points on the edges
//         I               -- current index
// 0 1 2 3 4 5 6 7 8 9 A   -- point index
// p p p p p p p           -- point
//   1 1 1 1               -- first derivative
//     2 2 2               -- second derivative
//     i i i               -- i curvature value
//       k k               -- smoothed curvature vale
//         c               -- class
// 0 1 2 3 4 5 6 7 8 9 A   -- point index
// ------------------------------------------------------------------------
TracePointClass InputSmoothing::calculate_class( size_t index ) {
    assert( can_be_calculated_now( index ) );
    if ( index == 0 ) {
        return TracePointClass::SharpEdge;
    }
    Point2 p2 = fetch_point( index - 2 );
    Point2 p3 = fetch_point( index - 1 );
    Point2 p4 = fetch_point( index );
    Point2 p5 = fetch_point( index + 1 );
    Point2 p6 = fetch_point( index + 2 );

    double curvature = calculate_i_curvature( p2, p3, p4, p5, p6 );
    if ( std::isnan( m_prev_curvature ) ) {
        m_prev_curvature = curvature;
    }
    auto prev_curvature = m_prev_curvature;
    m_prev_curvature = curvature;
    m_curvature[index] = curvature;

    if ( sign( prev_curvature ) != sign( curvature ) ) {
        return TracePointClass::Inflection;
    } else if ( curvature > 2.0 ) {
        return TracePointClass::SharpEdge;
    } else {
        return TracePointClass::Normal;
    }
}

void InputSmoothing::append_trace( const Trace &trace ) {
    for ( auto t : trace.samples ) {
        m_points_buffer.push_back( t.point );
    }
    m_curvature.resize( m_points_buffer.size() );
}

template <class T> void print_array( const vector<T> &v ) {
    for ( auto it = begin( v ); it != end( v ); ++it ) {
        printf( "%s%s", to_string( *it ).c_str(),
                it != end( v ) ? ", " : "\n" );
    }
}
