#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

#define EPS 0.000000001
#define PI 3.1415926535897932384626433832795028841971

using namespace std;

enum IntersectionType {NONE, POINT, LINE, PLANE};
enum SLENumSolutions {ZERO, ONE, INFINITE};

long double RadToDeg(long double rad) {
	return 180 / PI * rad;
}

long double DegToRad(long double deg) {
	return PI / 180 * deg;
}

struct Intersection {
	IntersectionType type;
	void* object;
};

class Point {
	public:
		vector<long double> coordinates;

	Point(long double x, long double y) {
		coordinates.push_back(x);
		coordinates.push_back(y);
	}

	Point(long double x, long double y, long double z) {
		coordinates.push_back(x);
		coordinates.push_back(y);
		coordinates.push_back(z);
	}

	Point(vector<long double> coords) {
		coordinates = coords;
	}

	Point(const Point& p) {
		coordinates = p.coordinates;
	}

	int GetRank() {
		return coordinates.size();
	}

	long double GetX() {
		return coordinates[0];
	}

	long double GetY() {
		return coordinates[1];
	}

	long double GetZ() {
		return coordinates[2];
	}

	long double Get(int index) {
		return coordinates[index];
	}

	void SetX(long double x) {
		coordinates[0] = x;
	}

	void SetY(long double y) {
		coordinates[1] = y;
	}

	void SetZ(long double z) {
		coordinates[2] = z;
	}

	void Set(int index, long double value) {
		coordinates[index] = value;
	}

	long double SqrMagnitude() {
		long double value = 0;
		for(int i = 0; i < GetRank(); i++) {
			value += coordinates[i] * coordinates[i];
		}
		return value;
	}

	long double Magnitude() {
		return sqrt(SqrMagnitude());
	}

	void Add(Point p) {
		for(int i = 0; i < GetRank(); i++) {
			coordinates[i] += p.coordinates[i];
		}
	}

	void Subtract(Point p) {
		for(int i = 0; i < GetRank(); i++) {
			coordinates[i] -= p.coordinates[i];
		}
	}

	void Scale(long double factor) {
		for(int i = 0; i < GetRank(); i++) {
			coordinates[i] *= factor;
		}
	}

	void Normalize() {
		if(abs(Magnitude()) < EPS) {
			return;
		}
		Scale(1 / Magnitude());
	}

	void Rotate2D(long double angle, Point center = Point(0, 0)) {
		Subtract(center);
		long double newX = GetX() * cos(angle) - GetY() * sin(angle);
		long double newY = GetX() * sin(angle) + GetY() * cos(angle);
		SetX(newX);
		SetY(newY);
		Add(center);
	}

	void Rotate3D(long double angle, Point normal, Point center = Point(0, 0, 0)) {
		Point unitNormal(normal);
		unitNormal.Normalize();
		Subtract(center);
		Point rotation = *this * cos(angle) + unitNormal.CrossProduct(*this) * sin(angle) + unitNormal * (unitNormal.DotProduct(*this)) * (1 - cos(angle));
		SetX(rotation.GetX());
		SetY(rotation.GetY());
		SetZ(rotation.GetZ());
		Add(center);
	}

	long double DotProduct(Point p) {
		long double value = 0;
		for(int i = 0; i < GetRank(); i++) {
			value += coordinates[i] * p.coordinates[i];
		}
		return value;
	}

	Point CrossProduct(Point p) {
		Point a(*this);
		Point b(p);
		Point result = Point(0, 0, 0);
		while(a.GetRank() < 3) {
			a.coordinates.push_back(0);
		}
		while(b.GetRank() < 3) {
			b.coordinates.push_back(0);
		}
		result.SetX(a.GetY() * b.GetZ() - a.GetZ() * b.GetY());
		result.SetY(a.GetZ() * b.GetX() - a.GetX() * b.GetZ());
		result.SetZ(a.GetX() * b.GetY() - a.GetY() * b.GetX());
		return result;
	}

	Point operator + (Point p) const {
		Point addition(coordinates);
		addition.Add(p);
		return addition;
	}

	Point operator - (Point p) const {
		Point subtraction(coordinates);
		subtraction.Subtract(p);
		return subtraction;
	}

	Point operator * (long double amount) const {
		Point multiplication(coordinates);
		multiplication.Scale(amount);
		return multiplication;
	}

	Point operator / (long double amount) const {
		Point division(coordinates);
		division.Scale(1 / amount);
		return division;
	}

	bool operator == (Point p) const {
		if(coordinates.size() != p.GetRank()) {
			return false;
		}
		for(int i = 0; i < coordinates.size(); i++) {
			if(abs(coordinates[i] - p.coordinates[i]) > EPS) {
				return false;
			}
		}
		return true;
	}

	bool operator < (Point p) const {
		for(int i = 0; i < coordinates.size(); i++) {
			if(abs(coordinates[i] - p.coordinates[i]) > EPS) {
				return coordinates[i] < p.coordinates[i];
			}
		}
		return false;
	}
};

struct SLESolution {
	SLENumSolutions numSolutions;
	vector<long double> solutionVector;
};

Point Zero(int rank) {
	vector<long double> coords;
	for(int i = 0; i < rank; i++) {
		coords.push_back(0);
	}
	return Point(coords);
}

bool CCW(Point a, Point b, Point c) {
	Point ba = a - b;
	Point bc = c - b;
	return ba.CrossProduct(bc).GetZ() > -EPS;
}

long double Angle(Point a, Point b, Point c) {
	Point ba = a - b;
	Point bc = c - b;
	return acos(ba.DotProduct(bc) / (ba.Magnitude() * bc.Magnitude()));
}

class Matrix {
	public:
		vector<vector<long double>> matrix;

	Matrix(int rows, int cols) {
		for(int y = 0; y < rows; y++) {
			vector<long double> submatrix;
			for(int x = 0; x < cols; x++) {
				submatrix.push_back(0);
			}
			matrix.push_back(submatrix);
		}
	}

	Matrix(vector<vector<long double>> m) {
		matrix = m;
	}

	Matrix(const Matrix& m) {
		matrix = m.matrix;
	}

	long double Get(int row, int col) {
		return matrix[row][col];
	}

	void Set(int row, int col, long double value) {
		matrix[row][col] = value;
	}

	void SwapRow(int a, int b) {
		for(int x = 0; x < GetColumns(); x++) {
			long double tmp = matrix[a][x];
			matrix[a][x] = matrix[b][x];
			matrix[b][x] = tmp;
		}
	}

	void SwapColumn(int a, int b) {
		for(int y = 0; y < GetRows(); y++) {
			long double tmp = matrix[y][a];
			matrix[y][a] = matrix[y][b];
			matrix[y][b] = tmp;
		}
	}

	void DeleteRow(int row) {
		matrix.erase(matrix.begin() + row);
	}

	void DeleteColumn(int col) {
		for(int y = 0; y < GetRows(); y++) {
			matrix[y].erase(matrix[y].begin() + col);
		}
	}

	void InsertRow(int row) {
		vector<long double> newRow;
		for(int x = 0; x < GetColumns(); x++) {
			newRow.push_back(0);
		}
		matrix.insert(matrix.begin() + row, newRow);
	}

	void InsertColumn(int col) {
		for(int y = 0; y < GetRows(); y++) {
			matrix[y].insert(matrix[y].begin() + col, 0);
		}
	}

	int GetRows() {
		return matrix.size();
	}

	int GetColumns() {
		if(matrix.size() > 0) {
			return matrix[0].size();
		}
		return 0;
	}

	Matrix GetTranspose() {
		Matrix transpose(GetColumns(), GetRows());
		for(int y = 0; y < GetRows(); y++) {
			for(int x = 0; x < GetColumns(); x++) {
				transpose.Set(x, y, Get(y, x));
			}
		}
		return transpose;
	}

	Matrix GetSubmatrix(int rows, int cols, int startRow = 0, int startCol = 0) {
		Matrix submatrix(rows, cols);
		for(int y = 0; y < rows; y++) {
			for(int x = 0; x < cols; x++) {
				submatrix.matrix[y][x] = matrix[y + startRow][x + startCol];
			}
		}
		return submatrix;
	}

	Point Multiply(Point Xvec) {
		if(GetColumns() != Xvec.GetRank()) {
			cerr << "Invalid matrix/vector size for transformation" << endl;
			return Xvec;
		} else {
			vector<long double> components;
			for(int y = 0; y < GetRows(); y++) {
				long double comp = 0;
				for(int x = 0; x < GetColumns(); x++) {
					comp += Xvec.coordinates[x] * matrix[y][x];
				}
				components.push_back(comp);
			}
			return Point(components);
		}
	}

	long double GetDeterminant() {
		long double determinant = 0;
		if(GetRows() != GetColumns()) {
			cerr << "Cannot obtain the determinant of a non-square matrix" << endl;
			return 0;
		}
		if(GetRows() > 1) {
			for(int x = 0; x < GetColumns(); x++) {
				Matrix submatrix = GetSubmatrix(GetRows() - 1, GetColumns(), 1);
				submatrix.DeleteColumn(x);
				if(x % 2 == 0) {
					determinant += matrix[0][x] * submatrix.GetDeterminant();
				} else {
					determinant -= matrix[0][x] * submatrix.GetDeterminant();
				}
			}
		} else {
			return matrix[0][0];
		}
		return determinant;
	}

	SLESolution SolveEquations(Point Yvec) {
		SLESolution none;
		none.numSolutions = ZERO;
		SLESolution infinite;
		infinite.numSolutions = INFINITE; 
		if(Yvec.GetRank() != GetRows()) {
			return none;
		}
		Matrix reducedSystem(*this);
		reducedSystem.InsertColumn(reducedSystem.GetColumns());
		for(int y = 0; y < Yvec.GetRank(); y++) {
			reducedSystem.Set(y, reducedSystem.GetColumns() - 1, Yvec.coordinates[y]);
		}
		for(int y = reducedSystem.GetRows() - 1; y >= 0; y--) {
			if(Point(reducedSystem.matrix[y]) == Zero(reducedSystem.GetColumns())) {
				reducedSystem.DeleteRow(y);
			}
		}
		for(int a = reducedSystem.GetRows() - 1; a >= 0; a--) {
			for(int b = a - 1; b >= 0; b--) {
				Point pa(reducedSystem.matrix[a]);
				Point pb(reducedSystem.matrix[b]);
				pa.Normalize();
				pb.Normalize();
				if((pa == pb) || (pa == (pb * -1))) {
					reducedSystem.DeleteRow(a);
					break;
				}
			}
		}
		for(int a = reducedSystem.GetRows() - 1; a >= 0; a--) {
			for(int b = a - 1; b >= 0; b--) {
				Point pa(reducedSystem.matrix[a]);
				Point pb(reducedSystem.matrix[b]);
				pa.coordinates.pop_back();
				pb.coordinates.pop_back();
				pa.Normalize();
				pb.Normalize();
				if((pa == pb) || (pa == (pb * -1))) {
					return none;
				}
			}
		}
		if(reducedSystem.GetRows() < reducedSystem.GetColumns() - 1) {
			return infinite;
		}
		vector<long double> newYVec;
		vector<Point> XVecs;
		for(int y = 0; y < reducedSystem.GetRows(); y++) {
			newYVec.push_back(reducedSystem.Get(y, reducedSystem.GetColumns() - 1));
		}
		reducedSystem.DeleteColumn(reducedSystem.GetColumns() - 1);
		for(int k = 0; k < reducedSystem.GetRows() - 1; k++) {
			Matrix subsystem = reducedSystem.GetSubmatrix(reducedSystem.GetColumns(), reducedSystem.GetColumns(), k);
			vector<long double> subsystemYVecCoords;
			for(int i = 0; i < reducedSystem.GetColumns(); i++) {
				subsystemYVecCoords.push_back(newYVec[k + i]);
			}
			Point subsystemYVec(subsystemYVecCoords);
			Matrix LMatrix(subsystem.GetRows(), subsystem.GetColumns());
			for(int i = 0; i < subsystem.GetRows(); i++) {
				LMatrix.matrix[i][i] = 1;
			}
			Matrix UMatrix(subsystem);
			int swappingRow = 0;
			for(int x = 0; x < subsystem.GetColumns(); x++) {
				for(int y = swappingRow; y < subsystem.GetRows(); y++) {
					if(abs(UMatrix.matrix[y][x]) > EPS) {
						if(y != swappingRow) {
							UMatrix.SwapRow(y, swappingRow);
							long double tmp = subsystemYVec.coordinates[y];
							subsystemYVec.coordinates[y] = subsystemYVec.coordinates[swappingRow];
							subsystemYVec.coordinates[swappingRow] = tmp;
						}
						swappingRow++;
						break;
					}
				}
			}
			int currCol = 0;
			for(int y = 0; y < subsystem.GetRows(); y++) {
				int currRow = y;
				if(currCol >= subsystem.GetColumns()) {
					break;
				}
				while(UMatrix.matrix[currRow][currCol] == 0) {
					currRow++;
					if(currRow >= subsystem.GetRows()) {
						currRow = y;
						currCol++;
						if (currCol >= subsystem.GetColumns()) {
							break;
						}
					}
				}
				if(currCol < subsystem.GetColumns()) {
					if(currRow != y) {
						UMatrix.SwapRow(currRow, y);
					}
					for(int z = y + 1; z < subsystem.GetRows(); z++) {
						if(UMatrix.matrix[z][currCol] != 0) {
							long double mult = UMatrix.matrix[z][currCol] / UMatrix.matrix[y][currCol];
							LMatrix.matrix[z][currCol] = mult;
							for(int x = currCol; x < subsystem.GetColumns(); x++) {
								UMatrix.matrix[z][x] -= mult * UMatrix.matrix[y][x];
							}
						}
					}
					currCol++;
				}
			}
			vector<long double> Zvec;
			vector<long double> output;
			for(int y = 0; y < subsystem.GetRows(); y++) {
				long double current = subsystemYVec.coordinates[y];
				for(int x = 0; x < subsystem.GetColumns(); x++) {
					if(x < Zvec.size()) {
						current -= LMatrix.matrix[y][x] * Zvec[x];
					}
				}
				Zvec.push_back(current);
			}
			for(int y = subsystem.GetRows() - 1; y >= 0; y--) {
				long double current = Zvec[y];
				for(int x = subsystem.GetColumns() - 1; x >= 0; x--) {
					if(x > y) {
						current -= output[subsystem.GetColumns() - x - 1] * UMatrix.matrix[y][x];
					}
					if(x == y) {
						current /= UMatrix.matrix[y][x];
					}
				}
				output.push_back(current);
			}
			vector<long double> reverseOutput;
			while(output.size() > 0) {
				reverseOutput.push_back(output.back());
				output.pop_back();
			}
			Point XVec(reverseOutput);
			XVecs.push_back(XVec);
		}
		bool valid = true;
		for(int y = 1; y < XVecs.size(); y++) {
			if(!(XVecs[y] == XVecs[y - 1])) {
				valid = false;
			}
		}
		if(valid) {
			SLESolution one;
			one.numSolutions = ONE;
			one.solutionVector = XVecs[0].coordinates;
			return one;
		}
		return none;
	}
};

class Line {
	public:
		Point start;
		Point end;

	Line(Point s, Point e) : start(s), end(e) {}
	Line(long double sx, long double sy, long double ex, long double ey) : start(sx, sy), end(ex, ey) {}
	Line(long double sx, long double sy, long double sz, long double ex, long double ey, long double ez) : start(sx, sy, sz), end(ex, ey, ez) {}

	long double GetLength() {
		Point vec = end - start;
		return vec.Magnitude();
	}

	Point GetNearestPoint(Point p) {
		Point vec = end - start;
		Point lvec = p - start;
		long double n = vec.DotProduct(lvec) / (vec.Magnitude() * vec.Magnitude());
		return GetInterpolation(n);
	}

	Point GetNearestPointInSegment(Point p) {
		Point vec = end - start;
		Point lvec = p - start;
		long double n = vec.DotProduct(lvec) / (vec.Magnitude() * vec.Magnitude());
		if(n < 0) {
			n = 0;
		}
		if(n > 1) {
			n = 1;
		}
		return GetInterpolation(n);
	}

	Point GetVector() {
		return end - start;
	}

	bool IsCollinear(Point p) {
		return GetNearestPoint(p) == p;
	}

	bool IsInSegment(Point p) {
		if(IsCollinear(p)) {
			Point vec = end - start;
			Point lvec = p - start;
			long double n = vec.DotProduct(lvec) / (vec.Magnitude() * vec.Magnitude());
			if(n >= 0 && n <= 1) {
				return true;
			}
		}
		return false;
	}

	bool IsParallel(Line l) {
		Point vec = end - start;
		Point lvec = l.end - l.start;
		vec.Normalize();
		lvec.Normalize();
		return (vec == lvec) || (vec == (lvec * -1));
	}

	bool IsPerpendicular(Line l) {
		Point vec = end - start;
		Point lvec = l.end - l.start;
		return abs(vec.DotProduct(lvec)) < EPS;
	}

	bool Equals(Line l) {
		return IsCollinear(l.start) && IsCollinear(l.end);
	}

	bool EqualsSegment(Line l) {
		return (start == l.start && end == l.end) || (start == l.end && end == l.start);
	}

	Intersection GetIntersection(Line l) {
		Intersection none;
		none.type = NONE;
		if(IsParallel(l)) {
			if(IsCollinear(l.start)) {
				Line* linePtr = new Line(start, end);
				Intersection line;
				line.type = LINE;
				line.object = (void*)linePtr;
				return line;
			} else {
				return none;
			}
		}
		Matrix coefficients(start.GetRank(), 2);
		for(int y = 0; y < start.GetRank(); y++) {
			coefficients.Set(y, 0, (end - start).coordinates[y]);
			coefficients.Set(y, 1, -(l.end - l.start).coordinates[y]);
		}
		SLESolution solution = coefficients.SolveEquations(l.start - start);
		if(solution.numSolutions == ONE) {
			long double s = solution.solutionVector[0];
			Point* pointPtr = new Point(GetInterpolation(s));
			Intersection point;
			point.type = POINT;
			point.object = (void*)pointPtr;
			return point;
		}
		return none;
	}

	Intersection GetSegmentIntersection(Line l) {
		Intersection none;
		none.type = NONE;
		Intersection intersect = GetIntersection(l);
		if(intersect.type == POINT) {
			Point p = *(Point*)intersect.object;
			if(IsInSegment(p) && l.IsInSegment(p)) {
				return intersect;
			}
		} else if (intersect.type == LINE) {
			vector<Point> pointsInSegment;
			vector<Point> uniquePoints;
			if(IsInSegment(l.start)) {
				pointsInSegment.push_back(l.start);
			}
			if(IsInSegment(l.end)) {
				pointsInSegment.push_back(l.end);
			}
			if(l.IsInSegment(start)) {
				pointsInSegment.push_back(start);
			}
			if(l.IsInSegment(end)) {
				pointsInSegment.push_back(end);
			}
			for(int i = 0; i < pointsInSegment.size(); i++) {
				bool unique = true;
				for(int j = 0; j < uniquePoints.size(); j++) {
					if(pointsInSegment[i] == uniquePoints[j]) {
						unique = false;
						break;
					}
				}
				if(unique) {
					uniquePoints.push_back(pointsInSegment[i]);
				}
			}
			if (uniquePoints.size() == 1) {
				Point* pointPtr = new Point(uniquePoints[0]);
				Intersection pointIntersect;
				pointIntersect.type = POINT;
				pointIntersect.object = (void*)pointPtr;
				return pointIntersect;
			} else if (uniquePoints.size() == 2) {
				Line* linePtr = new Line(uniquePoints[0], uniquePoints[1]);
				Intersection lineIntersect;
				lineIntersect.type = LINE;
				lineIntersect.object = (void*)linePtr;
				return lineIntersect;
			}
		}
		return none;
	}

	Point GetInterpolation(long double fraction) {
		Point vec = end - start;
		vec.Scale(fraction);
		vec.Add(start);
		return vec;
	}
};

class Plane {
	public:
		Point center;
		Point normal;

	Plane(Point c, Point a, Point b) : center(c), normal((a - c).CrossProduct(b - c)) {}
	Plane(Point c, Point n) : center(c), normal(n) {}

	bool ContainsPoint(Point p) {
		return abs(normal.DotProduct(p - center)) <= EPS;
	}

	bool ContainsLine(Line l) {
		return ContainsPoint(l.start) && ContainsPoint(l.end);
	}

	bool IsParallel(Line l) {
		return abs(normal.DotProduct(l.end - l.start)) <= EPS;
	}

	bool IsParallel(Plane p) {
		Point unitNormal(normal);
		Point otherNormal(p.normal);
		unitNormal.Normalize();
		otherNormal.Normalize();
		return (unitNormal == otherNormal) || (unitNormal == (otherNormal * -1));
	}

	bool IsPerpendicular(Line l) {
		Point unitNormal(normal);
		Point lineNormal = l.end - l.start;
		unitNormal.Normalize();
		lineNormal.Normalize();
		return (unitNormal == lineNormal) || (unitNormal == (lineNormal * -1));
	}

	bool IsPerpendicular(Plane p) {
		return abs(normal.DotProduct(p.normal)) <= EPS;
	}

	Intersection GetIntersection(Line l) {
		Intersection none;
		none.type = NONE;
		if(IsParallel(l)) {
			if(ContainsLine(l)) {
				Line* linePtr = new Line(l.start, l.end);
				Intersection lineIntersect;
				lineIntersect.type = LINE;
				lineIntersect.object = (void*)linePtr;
				return lineIntersect;
			} else {
				return none;
			}
		} else {
			long double distance = (center - l.start).DotProduct(normal) / (l.end - l.start).DotProduct(normal);
			Point* pointPtr = new Point(l.GetInterpolation(distance));
			Intersection pointIntersect;
			pointIntersect.type = POINT;
			pointIntersect.object = (void*)pointPtr;
			return pointIntersect;
		}
	}

	Intersection GetIntersection(Plane p) {
		Intersection none;
		none.type = NONE;
		if(IsParallel(p)) {
			if(ContainsPoint(p.center)) {
				Plane* planePtr = new Plane(center, normal);
				Intersection planeIntersect;
				planeIntersect.type = PLANE;
				planeIntersect.object = (void*)planePtr;
				return planeIntersect;
			} else {
				return none;
			}
		} else {
			Point lineVec = normal.CrossProduct(p.normal);
			Matrix coefficients(2, normal.GetRank());
			vector<long double> columnIndices;
			Point YVec(normal.DotProduct(center), p.normal.DotProduct(p.center));
			for(int x = 0; x < coefficients.GetColumns(); x++) {
				coefficients.Set(0, x, normal.Get(x));
				coefficients.Set(1, x, p.normal.Get(x));
				columnIndices.push_back(x);
			}
			for(int x = coefficients.GetColumns() - 1; x >= 0; x--) {
				if(coefficients.Get(0, x) == 0 && coefficients.Get(1, x) == 0) {
					coefficients.DeleteColumn(x);
					columnIndices.erase(columnIndices.begin() + x);
				}
			}
			for(int x = coefficients.GetColumns() - 1; x >= 0; x--) {
				if(coefficients.Get(0, x) == coefficients.Get(1, x) && coefficients.GetColumns() > 2) {
					coefficients.DeleteColumn(x);
					columnIndices.erase(columnIndices.begin() + x);
				}
			}
			while(coefficients.GetColumns() > 2) {
				coefficients.DeleteColumn(coefficients.GetColumns() - 1);
				columnIndices.pop_back();
			}
			SLESolution XVec = coefficients.SolveEquations(YVec);
			if(XVec.numSolutions == ONE) {
				Point lineStart = Zero(normal.GetRank());
				for(int x = 0; x < XVec.solutionVector.size(); x++) {
					lineStart.Set(columnIndices[x], XVec.solutionVector[x]);
				}
				Line* linePtr = new Line(lineStart, lineStart + lineVec);
				Intersection lineIntersect;
				lineIntersect.type = LINE;
				lineIntersect.object = (void*)linePtr;
				return lineIntersect;
			}
			return none;
		}
	}

	Point GetNearestPoint(Point p) {
		Intersection intersect = GetIntersection(Line(p, p + normal));
		if(intersect.type == POINT) {
			return *(Point*)intersect.object;
		} else {
			cerr << "Error finding nearest point to plane" << endl;
			return p;
		}
	}
};

class NSphere {
	public:
		Point center;
		long double radius;

	NSphere(Point c, long double r) : center(c) {
		radius = r;
	}

	NSphere(long double x, long double y, long double r) : center(x, y) {
		radius = r;
	}

	NSphere(long double x, long double y, long double z, long double r) : center(x, y, z) {
		radius = r;
	}

	bool ContainsPoint(Point p) {
		Point rp = p - center;
		return rp.SqrMagnitude() - radius * radius <= EPS;
	}

	bool PointOnBorder(Point p) {
		Point rp = p - center;
		return abs(rp.SqrMagnitude() - radius * radius) <= EPS;
	}

	bool Intersects(NSphere s) {
		Point rr = s.center - center;
		return rr.Magnitude() - radius - s.radius <= EPS;
	}

	bool FullyContains(NSphere s) {
		Point rr = s.center - center;
		return rr.Magnitude() + s.radius - radius <= EPS;
	}

	long double GetSurfaceArea() {
		long double rank = center.GetRank();
		return 2 * pow(PI, rank / 2.0) * pow(radius, rank - 1) / tgamma(rank / 2.0);
	}

	long double GetVolume() {
		long double rank = center.GetRank();
		return pow(PI, rank / 2.0) * pow(radius, rank) / tgamma(rank / 2.0 + 1);
	}

	Point GetNearestPointOnSurface(Point p) {
		Point vec = p - center;
		vec.Normalize();
		vec.Scale(radius);
		vec.Add(center);
		return vec;
	}
};

class Polygon {
	public:
		vector<Point> vertices;

	Polygon(vector<Point> points) {
		vertices = points;
	}

	long double GetArea() {
		long double area = 0;
		for(int i = 0; i < vertices.size(); i++) {
			area += vertices[i].GetX() * vertices[(i + 1) % vertices.size()].GetY() - vertices[(i + 1) % vertices.size()].GetX() * vertices[i].GetY();
		}
		return abs(area) / 2;
	}

	long double GetPerimeter() {
		long double perimeter = 0;
		for(int i = 0; i < vertices.size(); i++) {
			perimeter += GetEdge(i).GetLength();
		}
		return perimeter;
	}

	bool IsConvex() {
		RemoveDuplicates();
		bool directionSet = false;
		bool direction = false;
		for(int i = 0; i < vertices.size(); i++) {
			Point v1 = vertices[i];
			Point v2 = vertices[(i + 1) % vertices.size()];
			Point v3 = vertices[(i + 2) % vertices.size()];
			if(abs(Angle(v1, v2, v3) - PI) < EPS) {
				continue;
			}
			if(!directionSet) {
				directionSet = true;
				direction = CCW(v1, v2, v3);
			} else {
				if(CCW(v1, v2, v3) != direction) {
					return false;
				}
			}
		}
		return true;
	}

	bool ContainsPoint(Point p) {
		if(PointOnBorder(p)) {
			return true;
		}
		long double angleSum = 0;
		for(int i = 0; i < vertices.size(); i++) {
			Point v1 = vertices[i];
			Point v2 = vertices[(i + 1) % vertices.size()];
			if(CCW(v1, p, v2)) {
				angleSum += Angle(v1, p, v2);
			} else {
				angleSum -= Angle(v1, p, v2);
			}
		}
		return abs(abs(angleSum) - 2 * PI) < EPS;
	}

	bool PointOnBorder(Point p) {
		for(int i = 0; i < vertices.size(); i++) {
			if(GetEdge(i).IsInSegment(p)) {
				return true;
			}
		}
		return false;
	}

	void RemoveDuplicates() {
		for(int i = vertices.size() - 1; i >= 0; i--) {
			if(vertices[i] == vertices[(i + 1) % vertices.size()]) {
				vertices.erase(vertices.begin() + i);
			}
		}
	}

	void RemoveCollinearPoints() {
		for(int i = vertices.size() - 1; i >= 0; i--) {
			if(Line(vertices[(i - 1 + vertices.size()) % vertices.size()], vertices[(i + 1) % vertices.size()]).IsCollinear(vertices[i])) {
				vertices.erase(vertices.begin() + i);
			}
		}
	}

	Line GetEdge(int index) {
		Point v1(vertices[index % vertices.size()]);
		Point v2(vertices[(index + 1) % vertices.size()]);
		Line edge(v1, v2);
		return edge;
	}

	Polygon GetConvexHull() {
		int lowest = 0;
		for(int i = 1; i < vertices.size(); i++) {
			if(vertices[i].GetY() < vertices[lowest].GetY()) {
				lowest = i;
			} else if (vertices[i].GetY() == vertices[lowest].GetY()) {
				if(vertices[i].GetX() < vertices[lowest].GetX()) {
					lowest = i;
				}
			}
		}
		vector<Point> sortedVertices;
		for(int i = 0; i < vertices.size(); i++) {
			if(i != lowest) {
				sortedVertices.push_back(vertices[i]);
			}
		}
		Point lowestPoint = vertices[lowest];
		Point lowestDX = lowestPoint + Point(1, 0);
		sort(sortedVertices.begin(), sortedVertices.end(), [lowestPoint, lowestDX](Point a, Point b) -> bool {
			if(Angle(a, lowestPoint, lowestDX) < Angle(b, lowestPoint, lowestDX)) {
				return true;
			} else if (Angle(a, lowestPoint, lowestDX) == Angle(b, lowestPoint, lowestDX)) {
				return (a - lowestPoint).Magnitude() < (b - lowestPoint).Magnitude();
			}
			return false;
		});
		vector<Point> hull;
		hull.push_back(lowestPoint);
		hull.push_back(sortedVertices[0]);
		for(int i = 1; i < sortedVertices.size(); i++) {
			while(!CCW(sortedVertices[i], hull[hull.size() - 1], hull[hull.size() - 2])) {
				hull.pop_back();
			}
			hull.push_back(sortedVertices[i]);
		}
		Polygon convexHull(hull);
		return convexHull;
	}
};

class Triangle : public Polygon {
	public:
		Triangle(Point a, Point b, Point c) : Polygon({a, b, c}) {}

	bool HasRightAngle() {
		Line l1 = GetEdge(0);
		Line l2 = GetEdge(1);
		Line l3 = GetEdge(2);
		return l1.IsPerpendicular(l2) || l1.IsPerpendicular(l3) || l2.IsPerpendicular(l3);
	}

	int GetEqualSides() {
		Line l1 = GetEdge(0);
		Line l2 = GetEdge(1);
		Line l3 = GetEdge(2);
		int equalSides = 0;
		if(abs(l1.GetLength() - l2.GetLength()) < EPS) {
			equalSides++;
		};
		if(abs(l1.GetLength() - l3.GetLength()) < EPS) {
			equalSides++;
		};
		if(abs(l2.GetLength() - l3.GetLength()) < EPS) {
			equalSides++;
		};
		return equalSides;
	}

	NSphere GetCircumcircle() {
		Point midpointA = GetEdge(0).GetInterpolation(0.5);
		Point midpointB = GetEdge(1).GetInterpolation(0.5);
		Point perpA = GetEdge(0).GetVector();
		Point perpB = GetEdge(1).GetVector();
		perpA.Rotate2D(DegToRad(90));
		perpB.Rotate2D(DegToRad(90));
		Line bisectorA(midpointA, midpointA + perpA);
		Line bisectorB(midpointB, midpointB + perpB);
		Intersection intersection = bisectorA.GetIntersection(bisectorB);
		Point center = *(Point*)intersection.object;
		return NSphere(center, (center - vertices[0]).Magnitude());
	}

	NSphere GetIncircle() {
		long double a = GetEdge(1).GetLength();
		long double b = GetEdge(2).GetLength();
		long double c = GetEdge(0).GetLength();
		Point midpoint(a * vertices[0].GetX() + b * vertices[1].GetX() + c * vertices[2].GetX(), a * vertices[0].GetY() + b * vertices[1].GetY() + c * vertices[2].GetY());
		midpoint.Scale(1 / (a + b + c));
		long double radius = 2 * GetArea() / GetPerimeter();
		NSphere incircle(midpoint, radius);
		return incircle;
	}
};

class Quaternion {
	public:
		long double real;
		long double x;
		long double y;
		long double z;

	Quaternion(long double qr, long double qx, long double qy, long double qz) {
		real = qr;
		x = qx;
		y = qy;
		z = qz;
	}

	Quaternion(long double r, Point v) {
		real = r;
		x = v.GetX();
		y = v.GetY();
		z = v.GetZ();
	}

	Quaternion(const Quaternion &q) {
		real = q.real;
		x = q.x;
		y = q.y;
		z = q.z;
	}

	long double GetReal() {
		return real;
	}

	Point GetImaginary() {
		Point p(x, y, z);
		return p;
	}

	long double GetX() {
		return x;
	}

	long double GetY() {
		return y;
	}

	long double GetZ() {
		return z;
	}

	long double SqrMagnitude() {
		return real * real + x * x + y * y + z * z;
	}

	long double Magnitude() {
		return sqrt(SqrMagnitude());
	}

	void Normalize() {
		double mag = Magnitude();
		if(abs(mag) < EPS) {
			return;
		}
		real /= mag;
		x /= mag;
		y /= mag;
		z /= mag;
	}

	Quaternion Conjugate() {
		Quaternion conjugate(real, -x, -y, -z);
		return conjugate;
	}

	Quaternion Inverse() {
		long double sqr = SqrMagnitude();
		Quaternion inverse(real / sqr, -x / sqr, -y / sqr, -z / sqr);
		return inverse;
	}

	Quaternion operator + (Quaternion q) const {
		Quaternion addition(real + q.real, x + q.x, y + q.y, z + q.z);
		return addition;
	}

	Quaternion operator - (Quaternion q) const {
		Quaternion subtraction(real - q.real, x - q.x, y - q.y, z - q.z);
		return subtraction;
	}

	Quaternion operator * (Quaternion q) const {
		Quaternion multiplication(real * q.real - x * q.x - y * q.y - z * q.z, real * q.x + x * q.real + y * q.z - z * q.y, real * q.y - x * q.z + y * q.real + z * q.x, real * q.z + x * q.y - y * q.x + z * q.real);
		return multiplication;
	}

	Quaternion operator / (Quaternion q) const {
		Quaternion current(real, x, y, z);
		return current * q.Inverse();
	}

	bool operator == (Quaternion q) const {
		return !(abs(real - q.real) > EPS || abs(x - q.x) > EPS || abs(y - q.y) > EPS || abs(z - q.z) > EPS);
	}
};