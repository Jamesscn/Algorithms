#include <utility>
#include <vector>
#include <string>

using namespace std;

/**
 * Algorithm that aligns two strings A and B by inserting spaces and deleting characters
 * Longest Common Subsequence can be found with the following result:
 * needlemanWunsch(a, b, true, '-', 1, -INF, 0)
 * @param a The first string
 * @param b The second string
 * @param firstResult Makes the function stop after finding the first alignment
 * @param empty The character to place when a space is inserted
 * @param match The score that is used when a match is found
 * @param mismatch The score that is used when a mismatch is found
 * @param indel The score that is used when a character is deleted or a space is inserted
 * @return A vector of pairs of best alignments
 */
vector<pair<string, string>> needlemanWunsch(string a, string b, bool firstResult = false, char empty = '-', int match = 1, int mismatch = -1, int indel = -1) {
	vector<pair<string, string>> output;
	int scores[a.size() + 1][b.size() + 1];
	int paths[a.size() + 1][b.size() + 1];
	for(int y = 0; y <= a.size(); y++) {
		scores[y][0] = y * indel;
		paths[y][0] = 1;
	}
	for(int x = 0; x <= b.size(); x++) {
		scores[0][x] = x * indel;
		paths[0][x] = 2;
	}
	paths[0][0] = 0;
	for(int y = 1; y <= a.size(); y++) {
		for(int x = 1; x <= b.size(); x++) {
			int scoreTop = scores[y - 1][x] + indel;
			int scoreLeft = scores[y][x - 1] + indel;
			int scoreDiagonal = scores[y - 1][x - 1];
			if(a[y - 1] == b[x - 1]) {
				scoreDiagonal += match;
			} else {
				scoreDiagonal += mismatch;
			}
			scores[y][x] = max(max(scoreTop, scoreLeft), scoreDiagonal);
			paths[y][x] = 0;
			if(scoreTop == scores[y][x]) {
				paths[y][x] += 1;
			}
			if(scoreLeft == scores[y][x]) {
				paths[y][x] += 2;
			}
			if(scoreDiagonal == scores[y][x]) {
				paths[y][x] += 4;
			}
		}
	}
	pair<int, int> corner = make_pair(a.size(), b.size());
	vector<pair<int, int>> path;
	vector<pair<int, int>> movements;
	vector<int> distances;
	movements.push_back(corner);
	distances.push_back(0);
	while(movements.size() > 0) {
		pair<int, int> current = movements.back();
		int distance = distances.back();
		movements.pop_back();
		distances.pop_back();
		int y = current.first;
		int x = current.second;
		while(path.size() > distance) {
			path.pop_back();
		}
		path.push_back(current);
		int pathCode = paths[y][x];
		if(pathCode >= 4) {
			pathCode -= 4;
			movements.push_back(make_pair(y - 1, x - 1));
			distances.push_back(distance + 1);
		}
		if(pathCode >= 2) {
			pathCode -= 2;
			movements.push_back(make_pair(y, x - 1));
			distances.push_back(distance + 1);
		}
		if(pathCode >= 1) {
			movements.push_back(make_pair(y - 1, x));
			distances.push_back(distance + 1);
		}
		if(x == 0 && y == 0) {
			string newStringA = "";
			string newStringB = "";
			for(int i = path.size() - 2; i >= 0; i--) {
				int dy = path[i].first - path[i + 1].first;
				int dx = path[i].second - path[i + 1].second;
				if(dy == 1 && dx == 1) {
					newStringA += a[path[i].first - 1];
					newStringB += b[path[i].second - 1];
				} else if (dy == 1) {
					newStringA += a[path[i].first - 1];
					newStringB += empty;
				} else {
					newStringA += empty;
					newStringB += b[path[i].second - 1];
				}
			}
			output.push_back(make_pair(newStringA, newStringB));
			if(firstResult) {
				return output;
			}
		}
	}
	return output;
}