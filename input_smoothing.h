//
// Created by Lyubomyr Semkiv on 7/21/17.
//

#ifndef CURVES_INPUT_SMOOTHING_H
#define CURVES_INPUT_SMOOTHING_H

#include <iostream>
#include <vector>
#include <functional>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <sstream>
#include <iomanip>
#include <set>

#include "input_smoothing_types.h"

using namespace std;

// -----------------------------------------------------------------------

static const double DERIV_NAN = std::nan( "" );

// Expected that for each new trace approximation object is reset.
class InputSmoothing {
    static constexpr double default_test_tollerance = 10;

  public:
    InputSmoothing( double tollerance = default_test_tollerance);
    void on_trace_diff_available( Trace trace );
    void on_trace_end();
    void append_trace( const Trace &trace );

private:
    void process_next_point(size_t point_index, TracePointClass _class);

public:
    std::string fetch_curvature_around(size_t around, size_t back_n, size_t forth_n);

private:
    Point2 fetch_point(size_t at) const;
    size_t points_before_count(size_t index)const;
    size_t points_after_count(size_t index) const;
    bool can_be_calculated_later(size_t index) const;
    bool can_be_calculated_now(size_t index) const;
    TracePointClass calculate_class(size_t index) ;

  private:
    // classification data
    size_t m_offset = 0; // offset of fully processed points from buffer
    vector<Point2> m_points_buffer;
    double m_prev_curvature = 0.0;
    // simplification data
    Point2 m_corridor_p0;
    Point2 m_corridor_p1;
    size_t m_segment_begin_index = 0;
    size_t m_segment_end_index = 0;
    Point2 m_segment_tangent0; // begin of tangent line of current segment
    Point2 m_segment_tangent1; // end of tangent line of current segment
    // debug interface data
    vector<double> m_curvature;
};

#endif // CURVES_INPUT_SMOOTHING_H
