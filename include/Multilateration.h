// Łukasz Grochala
// http://people.duke.edu/~hpgavin/cee201/Nelder-Mead-2D.pdf
// http://www.jasoncantarella.com/downloads/NelderMeadProof.pdf

#ifndef MULTILATERATION_H
#define MULTILATERATION_H
#include <ArduinoSTL.h>
using namespace std;


namespace multilateration
{

const float epsilonX = 0.0000001; // epsilon for triangle size
const float epsilonF = 0.000001; // epsilon for error difference
const float epsilon = 0.0000000001;


template <typename T>
struct Point
{
    T x, y;
    Point(T x, T y);
    Point();
    Point<T> operator+ (const Point<T>& b);
    Point<T> operator- (const Point<T>& b);
    Point<T> operator* (const T& sc);
    Point<T> operator/ (const T& sc);
};

template <typename T>
Point<T> findSourceLocation (const std::vector<T>& distancesToStations, const std::vector<Point<T>>& stationsCoordinates);

template <typename T>
Point<T> findSourceLocation (const std::vector<T>& distancesToStations, const std::vector<Point<T>>& stationsCoordinates, Point<T> x0);
    
template <typename T>
Point<T> circularClipping (Point<T> a, const Point<T>& center, T radius);

template <typename T>
Point<T>::Point(T x, T y): x(x), y(y) {}

template <typename T>
Point<T>::Point() {}

template <typename T>
Point<T> Point<T>::operator+ (const Point<T>& b)
{
    return Point<T> (this->x + b.x, this->y + b.y);
}

template <typename T>
Point<T>  Point<T>::operator- (const Point<T>& b)
{
    return Point<T> (this->x - b.x, this->y - b.y);
}

template <typename T>
Point<T> Point<T>::operator* (const T& sc)
{
    return Point<T> (this->x * sc, this->y * sc);
}

template <typename T>
Point<T> Point<T>::operator/ (const T& sc)
{
    return Point<T> (this->x / sc, this->y / sc);
}

// Own implementation due to massive speedup in comparision to std::pow(x, 2) https://baptiste-wicht.com/posts/2017/09/cpp11-performance-tip-when-to-use-std-pow.html
template <typename T>
inline T pow2(T x)
{
    return (x*x);
}

template <typename T>
inline T _error(const Point<T>& a, const std::vector<T>& distancesToStations, const std::vector<Point<T>>& stationsCoordinates)
{
    T res = 0;
    auto itCoor = stationsCoordinates.begin();
    for (auto itDist = distancesToStations.begin();
        itDist != distancesToStations.end(); ++itDist, ++itCoor)
    {
        Point<T> dist (a.x - itCoor->x, a.y - itCoor->y);
        res += pow2(sqrt(pow2(dist.x) + pow2(dist.y))-*itDist);
    }
    return res;
}

template <typename T>
std::vector<Point<T>> _createTriangleFromInitialGuess (Point<T> x0)
{
    std::vector<Point<T>> triangle;
    triangle.reserve(3);
    triangle.push_back(x0);
    // creating 2 different points slightly moved in one dimension
    Point<T> x1 = x0;
    x1.x *= 1.05;
    if (x1.x == 0)
        x1.x = 0.00025;
    Point<T> x2 = x0;
    x2.y *= 1.05;
    if (x2.y == 0)
        x2.y = 0.00025;
    triangle.push_back(x1);
    triangle.push_back(x2);
    return triangle;
}


template <typename T>
std::vector <Point<T>> _findTriangleContainingEverySolution (const std::vector<T>& distancesToStations, const std::vector<Point<T>>& stationsCoordinates)
{
    // searching for rectangle that contains each possible solution
    T leftBound, topBound, bottomBound, rightBound;
    // finding rectangle that contains all stations
    leftBound = stationsCoordinates[0].x;
    rightBound = stationsCoordinates[0].x;
    topBound = stationsCoordinates[0].y;
    bottomBound = stationsCoordinates[0].y;
    for (const auto& station: stationsCoordinates)
    {
        // Min function is not used because of speed purposes - not sure about arduino compilator optimization
        if (station.x < leftBound)
            leftBound = station.x;
        if (station.x > rightBound)
            rightBound = station.x;
        if (station.y > topBound)
            topBound = station.y;
        if (station.y < bottomBound)
            bottomBound = station.y;
    }
    // finding maximum distance to station
    T maxDistance = *(max_element(distancesToStations.begin(), distancesToStations.end()));
    // extending rectangle to cover all possible source locations
    leftBound -= maxDistance;
    rightBound += maxDistance;
    topBound += maxDistance;
    bottomBound -= maxDistance;

    // finding triangle circumscribed about rectangle
    Point<T> x1, x2, x3;
    T rectangleWidth = rightBound - leftBound;
    T rectangleHeight = topBound - bottomBound;
    T averageSideLength = (rectangleHeight + rectangleWidth) / 2;
    x1.x = leftBound - averageSideLength;
    x1.y = bottomBound;
    x2.x = rightBound + averageSideLength;
    x2.y = bottomBound;
    x3.x = (rightBound + leftBound) / 2;
    // counted from similar triangle of known both side lengths TODO Ładny opis
    x3.y = topBound + ((rectangleWidth / 2 * rectangleHeight) / averageSideLength);
    return {x1, x2, x3};
}

// checking if best and worst node error is small enough to stop iteration
template <typename T>
inline bool _areFunctionsValuesSeperatedEnough(T bestError, T worstError)
{
    return abs((bestError - worstError) / (bestError + epsilon)) > epsilonF;
}

// checking if triangle area is small enough to stop iteration
template <typename T>
inline bool _isTriangleBigEnough(const Point<T>& A, const Point<T>& B, const Point<T>& C)
{
    return abs((A.x - C.x) * (B.y - A.y) - (A.x - B.x) * (C.y - A.y)) > epsilonX;
}

template <typename T>
inline void _sortNodesErrors(Point<T>& best, Point<T>& good, Point<T>& worst, T& bestError, T& goodError, T& worstError)
{
    if (goodError <= worstError && goodError <= bestError)
        {
            swap (goodError, bestError);
            swap (good.x, best.x);
            swap (good.y, best.y);
        }
        else if (worstError <= goodError && worstError <= bestError)
        {
            swap (worstError, bestError);
            swap (worst.x, best.x);
            swap (worst.y, best.y);
        }
        if (worstError < goodError)
        {
            swap (worstError, goodError);
            swap (worst.x, good.x);
            swap (worst.y, good.y);
        }
}

template <typename T>
Point<T> findSourceLocation (const std::vector<T>& distancesToStations, const std::vector<Point<T>>& stationsCoordinates)
{
    return _findSourceLocation (distancesToStations, stationsCoordinates, _findTriangleContainingEverySolution(distancesToStations, stationsCoordinates));
}

template <typename T>
Point<T> findSourceLocation (const std::vector<T>& distancesToStations, const std::vector<Point<T>>& stationsCoordinates, Point<T> x0)
{
    return _findSourceLocation (distancesToStations, stationsCoordinates, _createTriangleFromInitialGuess(x0));
}

template <typename T>
Point<T> _findSourceLocation (const std::vector<T>& distancesToStations, const std::vector<Point<T>>& stationsCoordinates, std::vector<Point<T>>&& triangle)
{
    // finding triangle that contains every possible solution
    Point<T> best = triangle[0];
    Point<T> good = triangle[1];
    Point<T> worst = triangle[2];
    T bestError = _error(best, distancesToStations, stationsCoordinates);
    T goodError = _error(good, distancesToStations, stationsCoordinates);
    T worstError = _error(worst, distancesToStations, stationsCoordinates);
    while (_isTriangleBigEnough(best, good, worst) || _areFunctionsValuesSeperatedEnough(bestError, worstError))
    {
        // sorting results for each node of the triangle
        _sortNodesErrors(best, good, worst, bestError, goodError, worstError);

        // creating Point in the middle of best and good
        Point<T> midPoint = (best + good) / 2;

        // creating point r as a reflection of worst
        Point<T> r = midPoint * 2 - worst;
        T errorR = _error(r, distancesToStations, stationsCoordinates);

        // checking if error for R is between error of worst and best
        if (bestError < errorR && errorR < goodError)
            worst = r;

        // extending R when error for R is the smallest one
        else if (errorR < bestError)
        {
            Point<T> e = r * 2 - midPoint;
            T errorE = _error(e, distancesToStations, stationsCoordinates);
            worst = errorE < errorR ? e : r;
        }
        else
        {
            // creating two points in 1/4 and 3/4 between worst and R and checking if some of them has smaller error than good
            Point<T> c1 = (midPoint + worst) / 2;
            Point<T> c2 = (midPoint + r) / 2;
            T c1Error = _error(c1, distancesToStations, stationsCoordinates);
            T c2Error = _error(c2, distancesToStations, stationsCoordinates);
            if (c2Error < c1Error)
            {
                c1 = c2;
                c1Error = c2Error;
            }
            if (c1Error < goodError)
                worst = c1;

            // shrinking triangle if none of the points had smaller error than good node
            else
            {
                good = (good + best) / 2;
                worst = (worst + best) / 2;
            }
        }
        bestError = _error(best, distancesToStations, stationsCoordinates);
        goodError = _error(good, distancesToStations, stationsCoordinates);
        worstError = _error(worst, distancesToStations, stationsCoordinates);
    }
    return best;
}

template <typename T>
Point<T> circularClipping (Point<T> a, const Point<T>& center, T radius)
{
    Point<T> clipped = a;
    // vector from Center to A point
    Point<T> centerAVector = a - center;
    T dist = sqrt(pow2(centerAVector.x) + pow2(centerAVector.y));
    if (dist > radius)
    {
        // counting vector from Center to closest point on a
        T ratio = radius / (dist);
        centerAVector = centerAVector * ratio;
        clipped = centerAVector + center;
    }
    return clipped;
}
    
}
    
#endif