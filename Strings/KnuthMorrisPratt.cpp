#include <string>
#include <vector>

using namespace std;

/**
 * Finds all occurrences of a pattern in a string of text
 * @param text The text to search
 * @param pattern The pattern to look for in the text
 * @return A vector containing all the indices of each occurrence
 */
vector<int> knuthMorrisPratt(string text, string pattern) {
	vector<int> output;
	int processTable[pattern.size()];
	for(int i = 0; i < pattern.size(); i++) {
		processTable[i] = 0;
	}
	int processLength = 0;
	int processIndex = 1;
	processTable[0] = 0;
	while(processIndex < pattern.size()) {
		if(pattern[processIndex] == pattern[processLength]) {
			processTable[processIndex] = processLength + 1;
			processLength++;
			processIndex++;
		} else {
			if(processLength == 0) {
				processTable[processIndex] = 0;
				processIndex++;
			} else {
				processLength = processTable[processLength - 1];
			}
		}
	}
	int textIndex = 0;
	int patternIndex = 0;
	while(textIndex < text.size()) {
		if(text[textIndex] == pattern[patternIndex]) {
			textIndex++;
			patternIndex++;
		} else {
			if(patternIndex == 0) {
				textIndex++;
			} else {
				patternIndex = processTable[patternIndex - 1];
			}
		}
		if(patternIndex == pattern.size()) {
			output.push_back(textIndex - patternIndex);
			patternIndex = processTable[patternIndex - 1];
		}
	}
	return output;
}