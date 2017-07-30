#pragma once

#include <vector>
#include <initializer_list>
#include <cstdint>

using std::vector;

struct Point2 {
    void set_invalid() { x = std::nan( "" ), y = std::nan( "" ); }
    bool is_invalid() const { return std::isnan( x ) && std::isnan( y ); }
    double x, y;
    bool operator==( Point2 other ) const {
        return ( x - x ) == 0.0 && ( y - y ) == 0.0;
    }
};


struct vec2d {
    double v[2];
    vec2d(double v0, double v1) { v[0] = v0; v[1] = v1; }
    vec2d(Point2 a, Point2 b) { v[0] = b.x - a.x; v[1] = b.y - a.y; }
    vec2d operator+(vec2d o) { return {v[0] + o.v[0], v[1] + o.v[1]}; }
    vec2d operator-(vec2d o) { return {v[0] - o.v[0], v[1] - o.v[1]}; }
    vec2d &operator+=(vec2d o) { v[0] += o.v[0], v[1] += o.v[1]; return *this; }
    vec2d &operator-=(vec2d o) { v[0] -= o.v[0], v[1] -= o.v[1]; return *this; }
    vec2d operator*(double s) const { return {v[0] * s, v[1] * s}; }
    vec2d operator/(double s) const { return {v[0] / s, v[1] / s}; }
    vec2d &operator*=(double s) { v[0] *= s, v[1] *= s; return *this; }
    vec2d &operator/=(double s) { v[0] /= s, v[1] /= s; return *this; }
    double operator[](size_t idx) const { return v[idx]; }
    double operator*(vec2d other) const { return v[0] * other.v[0] + v[1] * other.v[1]; }
};

namespace { double len(vec2d v) { return std::sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]); } }

static std::string to_string( Point2 p, bool long_format = false ) {
    std::stringstream s;
    s << std::setw(12);
    if (long_format)
        s << "Point2(" << p.x << ", " << p.y << ")";
    else
        s << p.x << ", " << p.y;
    return s.str();
}

struct Sample {
    Sample() : point(), timestamp( 0 ) {}
    Sample( Point2 point, uint64_t timestamp )
        : point( point ), timestamp( timestamp ) {}
    Point2 point;
    uint64_t timestamp;
};

struct Trace {
    Trace( std::initializer_list<Sample> l ) : samples( l ) {}
    Trace() {}
    vector<Sample> samples;
};

enum class PointClass { SharpEdge, Inflection };

struct InputSmootingResult {
    vector<Sample> trace;
    vector<Sample> simplified_trace;
    vector<double> first_derivatives;
    vector<double> second_derivatives;
    vector<double> curvature;
    vector<PointClass> classification;

    void check_invariant() {
        const size_t N = trace.size();
        assert( first_derivatives.size() == N );
        assert( second_derivatives.size() == N );
        assert( curvature.size() == N );
        assert( classification.size() == N );
    }
};

class InputSmoothingData {
  public:
    vector<Point2> simplification;
    vector<double> curvature;
    size_t iterative_offset;
    size_t iterative_count;
    size_t total_count;
};

using InputSmootingResultCb = std::function<void( InputSmootingResult )>;

// --------------------------------------------------------------

enum class TracePointClass { Normal, Inflection, SharpEdge };

class TraceAnalytics {
  public:
    double curvature_at( size_t index ) const {
        return _curvatures.at( index );
    }
    TracePointClass point_class_at( size_t index ) const {
        return _point_classes.at( index );
    }

    double derivative_at( int degree, size_t index ) const {
        return degree == 1 ? _first_derivatives.at( index )
                           : _second_derivatives.at( index );
    }

  private:
    vector<double> _first_derivatives;
    vector<double> _second_derivatives;
    vector<double> _curvatures;
    vector<TracePointClass> _point_classes;
};

enum class SampleClass { Regular, Simplified };

struct ExtendedSample {
    Point2 point;
    uint64_t timestamp;
    SampleClass kind;
    Sample as_sample() const { return Sample{point, timestamp}; }
};

struct SmoothingResult {
    SmoothingResult( TraceAnalytics analytics ) : analytics( analytics ) {}
    TraceAnalytics analytics;
    vector<Sample> original_trace;
    vector<Sample> simplified_trace;
};
