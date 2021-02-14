#include <utility>
#include <vector>
#include <deque>

using namespace std;

/*
Cross neighbour mask:
{
	make_pair(-1, 0),
	make_pair(0, 1),
	make_pair(1, 0),
	make_pair(0, -1)
}

Square neighbour mask:
{
	make_pair(-1, -1),
	make_pair(-1, 0),
	make_pair(-1, 1),
	make_pair(0, -1),
	make_pair(0, 1),
	make_pair(1, -1),
	make_pair(1, 0),
	make_pair(1, 1)
}
*/

/**
 * Function capable of doing both breadth first and depth first search on a 2D grid.
 * @param walls A two dimensional array of booleans containing walls the search cannot traverse.
 * @param rows The number of rows in the array.
 * @param cols The number of columns in the array.
 * @param start A pair of ints representing the starting Y and starting X of the search.
 * @param checkLocation A special function which is called at every step of the search, which takes in two pairs of ints,
 * the first one containing the current location and the second one containing the previous location.
 * If this function returns true, the search ends and the path to the final point is returned, otherwise an empty path is returned.
 * @param wrapAround If true, when the search reaches an edge it will wrap around to the other side, otherwise it will consider the edges to be walls.
 * @param breadthFirst If true, breadth first search is used, otherwise depth first search is applied.
 * @param neighbourMask A vector of relative points telling the search which adjacent cells can be considered neighbours.
 * @return A vector of points, which is empty if checkLocation returns false or the path to the end in reverse if it returns true.
 */
vector<pair<int, int>> Array2BDFS(bool** walls, int rows, int cols, pair<int, int> start, bool (*checkLocation)(pair<int, int>, pair<int, int>), bool wrapAround, bool breadthFirst, vector<pair<int, int>> neighbourMask) {
	vector<pair<int, int>> path;
	bool visited[rows][cols];
	pair<int, int> previous[rows][cols];
	for(int y = 0; y < rows; y++) {
		for(int x = 0; x < cols; x++) {
			visited[y][x] = walls[y][x];
			previous[y][x] = make_pair(y, x);
		}
	}
	deque<pair<int, int>> searchDeque;
	searchDeque.push_back(start);
	visited[start.first][start.second] = true;
	while(searchDeque.size() > 0) {
		pair<int, int> current;
		if(breadthFirst) {
			current = searchDeque.front();
			searchDeque.pop_front();
		} else {
			current = searchDeque.back();
			searchDeque.pop_back();
		}
		int y = current.first;
		int x = current.second;
		if(checkLocation(current, previous[y][x])) {
			while(current != start) {
				path.push_back(current);
				current = previous[current.first][current.second];
			}
			return path;
		}
		for(int i = 0; i < neighbourMask.size(); i++) {
			int newY = y + neighbourMask[i].first;
			int newX = x + neighbourMask[i].second;
			if(!wrapAround) {
				if(newY < 0 || newY >= rows) {
					continue;
				}
				if(newX < 0 || newX >= cols) {
					continue;
				}
			}
			newY = (newY % rows + rows) % rows;
			newX = (newX % cols + cols) % cols;
			if(!visited[newY][newX]) {
				previous[newY][newX] = current;
				visited[newY][newX] = true;
				searchDeque.push_back(make_pair(newY, newX));
			}
		}
	}
	return path;
}