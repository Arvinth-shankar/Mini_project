#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<iomanip>
#include<ctime>

constexpr double EPSILON {1e-4}; // tolerance for checking floating point numbers

// Class Point with point coordinates as data members
class Point{
    public : 
        double x_coord;
        double y_coord;

    public :

        Point() = default;
 
        Point(double x, double y) 
        : x_coord {x}, y_coord {y} {
        }

        bool operator==(const Point &p) const {
            return ((fabs(x_coord - p.x_coord) < EPSILON) && (fabs(y_coord - p.y_coord) < EPSILON));
        }

};

// defining struct Line with two points p1, p2
struct Line {
    Point p1;
    Point p2;
};

// holds the coefficients a,b,c of the line ax + by + c = 0 as data members
struct Line_equation {
    double coeff_a, coeff_b, coeff_c;
};

/* Function to find and return the central point for a given number of points in vector<Point>.
   Using Kahan summation algorithm to avoid errors in cummulative summation  */
Point find_center_point(const std::vector<Point> &points) {
    double centroid_x = 0.0, centroid_y = 0.0;
    double c_x = 0.0, c_y = 0.0;

    for (const auto &point : points) {
        double diff = point.x_coord - c_x;
        double t = centroid_x + diff;
        c_x = (t - centroid_x) - diff;
        centroid_x = t;

        diff = point.y_coord - c_y;
        t = centroid_y + diff;
        c_y = (t - centroid_y) - diff;
        centroid_y = t;
    }

    int n = points.size();
    return Point(centroid_x / n, centroid_y / n);
}

/* Reflect a point (p) over a line with two points p1 and p2.
   Result of this function will be a mirrored point over a line. */
Point mirrorpoint(const Point &p, const Line &line) {
    double xdiff = line.p2.x_coord - line.p1.x_coord; // vector x component of line
    double ydiff = line.p2.y_coord - line.p1.y_coord; // vector y component of line
    double xxp = p.x_coord - line.p1.x_coord;      
    double yyp = p.y_coord - line.p1.y_coord;
    double coeff = (xdiff * xxp + ydiff * yyp) / (xdiff * xdiff + ydiff * ydiff);
    double lx = line.p1.x_coord + xdiff * coeff; 
    double ly = line.p1.y_coord + ydiff * coeff; // lx and ly are the projected points
    return Point(2*lx-p.x_coord,2*ly-p.y_coord);
}

/* Function to check if a line reflects all points in the input point set 
   Returns true if all points are reflected over a line  */
bool does_line_reflect_all_points(const Line &line, const std::vector<Point> &points) {

    for (const auto &point : points) {
        Point reflected_point = mirrorpoint(point, line);
        if (std::find(points.begin(), points.end(), reflected_point) == points.end()) {
            return false;
        }
    }
    return true;
}

// Function to rotate a given point around the origin (0.0,0.0) by a specified angle in radians
Point rotate_point(const Point &point, double radians) {
    double cos_theta = cos(radians);
    double sin_theta = sin(radians);
    return (Point(point.x_coord * cos_theta - point.y_coord * sin_theta,
                  point.x_coord * sin_theta + point.y_coord * cos_theta));
}

/* Function to rotate an entire line around a specified center point by a given angle in radians
   Returns the rotated line indicated with two points p1, p2 */
Line rotate_line(const Line &line, double radians, const Point &center) {
    // Translate the line so that the center point becomes the origin    
    Point rotatedp1 = rotate_point(Point{line.p1.x_coord - center.x_coord, line.p1.y_coord - center.y_coord}, radians);
    Point rotatedp2 = rotate_point(Point{line.p2.x_coord - center.x_coord, line.p2.y_coord - center.y_coord}, radians);

    // Translate the rotated points back
    Point p1 {rotatedp1.x_coord + center.x_coord, rotatedp1.y_coord + center.y_coord};
    Point p2 {rotatedp2.x_coord + center.x_coord, rotatedp2.y_coord + center.y_coord};

    return Line{p1,p2};
}

/* Function to return a line perpendicular to line (specified by two points) and passing 
   through the midpoint of those two points  */
Line bisect_line(const Line &line) {
    Point midpoint = find_center_point({line.p1, line.p2});
    Point shifted_line = {line.p1.x_coord - midpoint.x_coord, line.p1.y_coord - midpoint.y_coord};
    Point shifted_bisector = {-shifted_line.y_coord, shifted_line.x_coord};

    return Line{midpoint, {shifted_bisector.x_coord + midpoint.x_coord, shifted_bisector.y_coord + midpoint.y_coord}};
}

/* Function to identify and return the first line that reflects all points
   Each pair of points has two possible lines of symmetry, one passing through those two points and
   the other perpendicular to the line and passing through midpoint of those two points */
Line identify_first_line(const std::vector<Point> &points, bool &first_line) {
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            Line passing_line = {points[i], points[j]};
            if (does_line_reflect_all_points(passing_line, points)) {
                first_line = true;
                return passing_line;
            }
            Line perpendicular_line = bisect_line(passing_line);
            if (does_line_reflect_all_points(perpendicular_line, points)) {
                first_line = true;
                return perpendicular_line;
            }
        }
    }
    first_line = false;
    return {{0, 0}, {0, 0}};
}

/* Function to remove duplicate coordinate points present in the input set (vector)
   This uses an std::vector and == operator in the Class Point to identify 
   the same coordinate points (using tolerance values) and remove them */
void remove_duplicates(std::vector<Point>& myVector) { 
    std::vector<Point> seen; 
  
    auto newEnd = remove_if( 
        myVector.begin(), myVector.end(), 
        [&seen](Point &point) { 
 
            if (std::find(seen.begin(), seen.end(), point) == seen.end()) { 
                seen.push_back(point); 
                return false;
            } 
            return true;
        }); 
  
    myVector.erase(newEnd, myVector.end()); 
} 

// Function to return all the factors for a given number 'n' 
std::vector<int> factorize(int n) {
    std::vector<int> result;
    for (int i = n; i >= 1; --i) {
        if (n % i == 0)
            result.push_back(i);
    }
    return result;
}

/* Function to identify all the possible lines of symmetry, if a symmetry exists for a given point set. 
   This function returns pair of points that each line of symmetry passes through.  */
auto find_lines(const std::vector<Point> &input_points) {

    const auto point_size = input_points.size();
    std::vector<Line> required_line {};

    // check the number of points in the set and find accordingly
    if(point_size ==0) {
        std::cout << "No points are present in the set!" << std::endl;
        return (required_line);
    }
    else if (point_size ==1) { // prints Infinite lines of symmetry in main() function
        return (required_line); 
    }
    else {
        const Point centroid = find_center_point(input_points);

        bool present_firstline = false;

        // Identify the first line of symmetry
        Line first_line = identify_first_line(input_points, present_firstline);

        // No symmetry exists if there is not a first line that reflects all points.
        if(!present_firstline) return (required_line);

        // Identify all the possible lines of symmetry by rotating the first line
        // Number of lines of symmetry is a factor of point_size or of point_size-1
        // Sort the factors in descending order
        std::vector<int> factors = factorize(point_size);
        std::vector<int> factors_N1 = factorize(point_size - 1);
        factors.insert(factors.end(), factors_N1.begin(), factors_N1.end());
        std::sort(factors.rbegin(), factors.rend());

        // Find the first highest factor, where all the points are reflected by the rotated line
        int factor = *std::find_if(factors.begin(), factors.end(), [&](int factor) {
            double radians = M_PI / factor;
            Line rotated_line = rotate_line(first_line, radians, centroid);
            return does_line_reflect_all_points(rotated_line, input_points);
            });
        
        // rotate the first line "factor" times to identify all the possible lines of symmetry
        for (int f = 0; f < factor; f++){
            required_line.push_back(rotate_line(first_line, f * M_PI / factor, centroid));
        }

    }

    return (required_line);
}

// Function to convert a line with two points to line with coefficients a, b, c (Line is represented as ax + by = c)
auto convert_line_points_to_coefficients(const Line &line){

    double coeff_a = line.p2.y_coord - line.p1.y_coord;
    double coeff_b = line.p1.x_coord - line.p2.x_coord;
    double coeff_c = coeff_a * (line.p1.x_coord) + coeff_b * (line.p1.y_coord);      
    
    if (fabs(coeff_a) < EPSILON) coeff_a = 0.0;
    if (fabs(coeff_b) < EPSILON) coeff_b = 0.0;
    if (fabs(coeff_c) < EPSILON) coeff_c = 0.0;

    return Line_equation{coeff_a, coeff_b, coeff_c};
}

int main() {
    
    time_t start, end;

    time(&start);

    int example_nr {};
    std::vector<Point> input_points;
    
    std::cout << "Enter an example to test: ";
    std::cin >> example_nr; 

    // Different test examples based on the user choice
    switch (example_nr)
    {
        case 1: // square (4 symmetry lines)
            input_points = {Point{-3.0,-6.0}, Point{-3.0,0.0}, Point {3.0,0.0}, 
            Point {3.0,-6.0}};
            break;
        case 2: // rectangle (2 symmetry lines)
            input_points = {Point{-2.0,3.0}, Point {-2.0,-3.0}, Point {2.0,3.0}, 
            Point {2.0,-3.0}};
            break;   
        case 3: // hexagon (6 symmetry lines)
            input_points = {Point{1.0,0.0}, Point{-1.0,0.0}, Point {0.5,0.8660}, 
            Point {-0.5,0.8660}, Point {0.5,-0.8660}, Point {-0.5,-0.8660}};
            break;
        case 4: // parallelogram (no symmetry lines)
            input_points = {Point{-2.0,2.0}, Point {-4.0,-2.0}, Point {2.0,-2.0}, 
            Point {4.0,2.0}};
            break;
        case 5: // parabola (1 symmetry line)
            input_points = {Point{-2.0,-1.0}, Point {-3.0,1.0}, Point {-4.0,3.0}, 
            Point {2.0,-1.0}, Point {3.0,1.0}, Point {4.0,3.0}};
            break;
        case 6: // equilateral triangle (3 symmetry lines)
            input_points = {Point{1.0,1.0}, Point {2.0,1.0}, Point {1.5,1.8660254037844386}};
            break;
        case 7: // just 6 points (2 symmetry lines)
            input_points = {Point{-2.0,3.0}, Point{-2.0,0.0}, Point {-2.0,-3.0}, 
            Point {2.0,3.0}, Point {2.0,0.0}, Point {2.0,-3.0}, Point {1.9999999,-3.00001}, Point{-2.0,0.0}};
            break; 
        case 8: // just 1 point (infinite symmetry lines)
            input_points = {Point{-2.0,3.0}};
            break;  
        case 9: // No points are present in the set
            input_points = {};
            break; 
        case 10: // 50000 points (2 lines of symmetry) - to check the time complexity of the current algorithm !
            // time taken for the whole program is also measured in main() function
            for (int i=0; i<50000; i++){
                input_points.push_back(Point{static_cast<double>(i), static_cast<double>(i)});
            }
            break;
        default:
            std::cout << "Enter a valid test case" << std::endl;
    }
    
    remove_duplicates(input_points);

    std::vector<Line> symmetry_lines = find_lines(input_points);

    /* If there is only one point in the input set, print
    infinite lines of symmetry, else print the number of symmetrical lines */
    if(input_points.size() == 1)
        std::cout << "Infinite lines of symmetry" << std::endl << std::endl;
    else 
        std::cout << "Number of lines of symmetry: " << symmetry_lines.size() << std::endl << std::endl;
    
    // Print the lines of symmetry in the format : ax + by = c
    if(symmetry_lines.size()) {
        std::cout << "The lines of symmetry are: " << std::endl << std::endl;

        Line_equation l;

        for (const auto &line : symmetry_lines){
            l = convert_line_points_to_coefficients(line);
            std::cout << std::setprecision(5);
            std::cout << l.coeff_a << "x + " << "(" << l.coeff_b << ")" << "y = " << l.coeff_c << std::endl << std::endl;
        }
    }
    
    time(&end);

    double time_taken = double(end - start); 
    std::cout << "Time taken by whole program is : " << std::fixed 
              << time_taken << std::setprecision(5); 
    std::cout << " sec " << std::endl; 
    return 0;
}
