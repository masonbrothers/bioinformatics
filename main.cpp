#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <sstream>
#include <cstdlib>

using namespace std;

struct Paladrome
{
	string sequence;
	unsigned int start;
	unsigned int stop;
};

struct Amino
{
	char peptide;
	string peptide3l;
	string peptideName;
	int number;
};

struct Codon
{
	string bases;
	char peptide;
	char type;
	bool special;
	bool start;
	bool stop;
	double bondAngle;
	double bondLength;
	string peptide3l;
	string peptideName;
};

struct Point
{
	double x;
	double y;
	double z;
	Codon type;
};

struct Frame
{
	vector<Codon> list;
	int startCodon;
	int lastCodon;
	string codingDNA;
	bool error;
	bool doesBreak;
	bool marker;
};

struct intext
{
	string code;
	bool extron;
};

struct SpliceSite
{
	//start[cut]first...last[cut]
	string start;
	string first;
	vector<string> middle;
	string last;
	string stop;
};

struct ThreeLettersCount
{
	string threeLetters;
	unsigned int number;
};

struct LetterCount
{
	char letter;
	unsigned int number;
};

struct Peptide
{
	Codon codon;
	vector<ThreeLettersCount> bases;
	unsigned int countBases()
	{
		unsigned int number = 0;
		for (int i = 0; i < bases.size(); i++)
		{
			number += bases.at(i).number;
		}
		return number;
	}

	vector<ThreeLettersCount> IUPACBases;
	unsigned int countIUPACBases()
	{
		unsigned int number = 0;
		for (int i = 0; i < bases.size(); i++)
		{
			number += IUPACBases.at(i).number;
		}
		return number;
	}
};

string RNAtoDNA(string RNA);

string DNAtoRNA(string DNA);

string DNAantiCodon(string DNA);

string RNAantiCodon(string RNA);

string flipString(string input);

vector<Codon> getCodonSheet(string filePath);

string getDNA(string filePath);

vector<string> parseStrings(string codonString);

Codon getAntiCodon(string inputCodon, vector<Codon> codonSheet);

Codon findCodon(const string &DNA, const vector<Codon> &codonSheet);

int findStartCodon(const string &DNA, const vector<Codon> &codonSheet, const int &start);

vector<Codon> getCodonSequence(int startCodon, int &lastCodon, const vector<Codon> &codonSheet, const string &DNA, string &codingDNA, bool &error, bool &doesBreak);

void printCodonSheetToConsole(const vector<Codon> &codonSheet);

bool printCodonSheetToFile(const vector<Codon> &codonSheet, const string &filePath);

vector<Frame> getFrames(const string &DNA, const vector<Codon> &codonSheet);

void printFramesToConsole(const vector<Frame> &Frames);

bool printFramesToFile(const vector<Frame> &Frames, const string &filePath);

bool printCodeToFile(string Code, string filePath);

void printCodeToConsole(string Code);

bool isPaladrome(const string Code);

vector<Paladrome> findPaladromes(const string &code, const int &minsize);

vector<Frame> findMinFrames(const vector<Frame> &Frames, const int &minsize);

vector<Frame> mergeFrames(const vector<Frame> &Frames1, const vector<Frame> &Frames2);

string DNAflipAntiString(string DNA);

bool printMachineFramesToFile(const vector<Frame> &Frames, const string &filePath);

bool isSubset(const char &large, const char &small);

double measureDistance(Point &point1, Point &point2);

double measureAngle(Point &point1, Point &point2, Point &point3);

Point translate(const Point &input, const double &x, const double &y, const double &z);

Point stretch(const Point &input, const double &x, const double &y, const double &z);

Point rotate_x(const Point &input, const double &angle);

Point rotate_y(const Point &input, const double &angle);

Point rotate_z(const Point &input, const double &angle);

vector<Amino> countPeptide(const vector<Codon> &codons);

vector<ThreeLettersCount> countCodons(const vector<Codon> &codons);

vector<ThreeLettersCount> countThreeLetters(const string &code, bool &error);

vector<LetterCount> countLetters(const string &code);

string getStatistics(const vector<Codon> &codonSheet, const string &DNA, bool &error, vector<Peptide> &peptides);

int getRandomNumber(const int &minimum, const int &maximum, const int &iteration);

ThreeLettersCount pickBase(vector<ThreeLettersCount> bases);

string getSequenceFromProtein(const Frame &proteinFrame, vector<Peptide> &peptides);

string getSequenceFromProtein(const string proteinString, vector<Peptide> &peptides);

char fromIUPAC(const char &large);

string fromIUPAC(const string &large);

Peptide getPeptideIUPACCodon(Peptide peptide);

vector<char> getMasters(const vector<char> &characters);

vector<Frame> getMinimumFrames(const string &DNA, const vector<Codon> &codonSheet, const unsigned int &minimumSize);

bool isDNA(const string &DNA);

bool isRNA(const string &RNA);

bool isNucleotide(const string &bases);

bool isPeptideLetter(const string &peptides, const vector<Codon> &codonSheet);

bool printToLog(string input);

string threeLetterPeptideToLetterPeptide(const string &peptides, vector<Codon> codonSheet);

string letterPeptideToThreeLetterPeptide(const string &peptides, vector<Codon> codonSheet);

bool isDNA(const string &DNA);

bool isRNA(const string &RNA);

bool isNucleotide(const string &bases);

bool isPeptideLetter(const string &peptides, const vector<Codon> &codonSheet);

bool isPeptideThreeLetters(const string &peptides, const vector<Codon> &codonSheet);

string getTypes(const string &code, const vector<Codon> &codonSheet);

string getTypes(const string &code, vector<int> &types, const vector<Codon> &codonSheet);

bool isLoop(const string &one, const string &two);

string findAndReplace(string input, const string &searchFor, const string &replaceWith, int &start);

string findAndReplaceOnceThrough(string input, const string &searchFor, const string &replaceWith);

string findAndReplace(string input, const string &searchFor, const string &replaceWith, int &start);

string findAndReplaceAll(string input, const string &searchFor, const string &replaceWith);


/*
This program does have some errors. Temporary fixes have been implemented. For some reason, the number 204 keeps appearing in the output. This sometimes manafests itself as the ascii 204 character.
This was fixed with the lines:

	if (tmpCodon.peptideName != "")
		{
			list.push_back(tmpCodon);
		}

Please fix before stable release.

Also a similar problem was found with the codon sheet. The last entry seems to be a partially complete entry of the last complete entry. This was solved by:
	if (peptideName != "")
	{
		codonSheet.push_back(tmp);
	}

*/

int main(int argc, char* argv[])
{
	vector<string> userConsoleArgs;
	for (int i = 0; i < argc; i++)
	{
		string temp(argv[i]);
		userConsoleArgs.push_back(temp);
	}

	return 1;
	bool helpFlag = false;
	bool errorFlag = false;
	for (int i = 0; i < userConsoleArgs.size(); i++)
	{
		if(userConsoleArgs.compare("-o") && (i + 1) != userConsoleArgs.size())
		{

		}
	}
	if (argc >= 1)
	bool error = 0;
	string codonSheetLocation = "CodonDatabase.txt";
	string DNAlocation = "DNAInput.txt";
	string outputLocation = "DNAOutput.txt";
	outputLocation = "CuratedOutput.txt";
	cout << "Copyright(C) 2016. Mason Brothers. All Rights Reserved." << endl;
	cout << "This program takes DNA and converts it to Proteins." << endl << endl;
	vector<Codon> codonSheet = getCodonSheet(codonSheetLocation);


	//printCodonSheetToConsole(codonSheet);

	string DNA = getDNA(DNAlocation);

	//string RNA = DNAtoRNA(DNA);
	//cout << RNA << endl;
	//printCodeToFile(DNA, outputLocation);
	//printCodeToConsole(DNA);
	//cout << "Error: " << error << endl << endl;

	int minProteinSize = 100;
	vector<Frame> bigFrames = getMinimumFrames(DNA, codonSheet, minProteinSize);
	printFramesToConsole(bigFrames);

	//printFramesToFile(bigFrames,"Frames.txt");

	cout << "There where " << bigFrames.size() << " reading frames with or more than " << minProteinSize << " amino acids.";
	cout << endl << endl;

	/*
	int start = 0;
	string test = "123Tes121233est123";
	string searchfor = "123";
	string replaceWith = "";
	//cout << findAndReplaceOnceThrough(test,searchfor,replaceWith) << endl;
	cout << findAndReplaceAll(test,searchfor,replaceWith) << endl;
	//cout << findAndReplace(test,searchfor,replaceWith,start) << endl;
	*/
	/*
	vector<Paladrome> paladromes;
	paladromes = findPaladromes(DNA, 2);
	for (int i = 0; i < paladromes.size(); i++)
	{
		cout << paladromes.at(i).sequence << " (" << paladromes.at(i).start << ", " << paladromes.at(i).stop << ")" << endl;
	}
	cout << endl;
	*/


	//printMachineFramesToFile(bigFrames, "Stuff.txt");

	//printFramesToConsole(bigFrames);


	//vector<Peptide> peptides;
	//cout << getStatistics(codonSheet, DNA, error, peptides);
	//cout << getTypes(DNA, codonSheet) << endl;


	/*
	vector<Frame> someFrames = findMinFrames(Frames,100);
	Frame test = someFrames.at(0);
	*/

	//string proteinString = "MKDLNNTKGNTKSEGSTERGNSGVDRGIVVPNTQIKMRFLNQVRYYSVNNNLKIGKDTNIELSKDTSTSDLLEFEKLVIDNINEENINNNLLSIIKNVDILILAYNRIKSKPGNITPGTTLETLDGINIIYLNKLSNELGTGKFKFKPMRIVNIPKPKGGIRPLSVGNPRDKIVQEVIRIILDTIFDKKISTHSHGFRKNISCQTAI";

	//string proteinString = "MV";
	//cout << proteinString << endl << getSequenceFromProtein(proteinString,peptides) << endl;

	/*
	forwardFrames = getFrames(sequence.str(),codonSheet);
	printFramesToConsole(forwardFrames);
	*/
	cout << endl;

	cout << "Error: " << error << endl << endl;

	/*
	int testNumber = 4;
	int testCount = 0;
	for (int i = 0; i < peptides.at(testNumber).bases.size(); i++)
	{
		cout << peptides.at(testNumber).bases.at(i).threeLetters << "\t" << peptides.at(testNumber).bases.at(i).number << endl;
		testCount += peptides.at(testNumber).bases.at(i).number;
	}
	cout << endl << "testCount: " << testCount << endl;
	cout << getRandomNumber(0, testCount, 1) << endl;
	cout << getRandomNumber(0, testCount, 1) << endl;
	cout << getRandomNumber(0, testCount, 1) << endl;
	cout << getRandomNumber(0, testCount, 1) << endl;
	cout << getRandomNumber(0, testCount, 1) << endl;
	cout << getRandomNumber(0, testCount, 1) << endl;
	cout << getRandomNumber(0, testCount, 1) << endl;
	cout << getRandomNumber(0, testCount, 1) << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	cout << pickBase(peptides.at(testNumber).bases).threeLetters << endl;
	*/
	/*
	Point one;
	one.x = 7;
	one.y = 7;
	one.z = 7;
	Point two;
	two.x = 1;
	two.y = 1;
	two.z = 1;
	Point three;
	three.x = -10;
	three.y = -10;
	three.z = -10;
	cout << measureAngle(one,two,three) << endl << measureDistance(one,two) << endl;
	*/
	//cout << rand() << "\t" << getRandomNumber(0,100,1) << endl;

	cout << "End of the Program!" << endl;
	return 0;
}

string motifToIUPAC(const string &input)
{
	//Note: This will only work for DNA not proteins
	string output = input;
	output = findAndReplaceAll(output,"[AT]","W");
	output = findAndReplaceAll(output,"[TA]","W");
	output = findAndReplaceAll(output,"[CG]","S");
	output = findAndReplaceAll(output,"[GC]","S");
	output = findAndReplaceAll(output,"[AC]","M");
	output = findAndReplaceAll(output,"[CA]","M");
	output = findAndReplaceAll(output,"[GT]","K");
	output = findAndReplaceAll(output,"[TG]","K");
	output = findAndReplaceAll(output,"[AG]","R");
	output = findAndReplaceAll(output,"[GA]","R");
	output = findAndReplaceAll(output,"[CT]","Y");
	output = findAndReplaceAll(output,"[TC]","Y");
	output = findAndReplaceAll(output,"[CGT]","B");
	output = findAndReplaceAll(output,"[CTG]","W");
	output = findAndReplaceAll(output,"[TCG]","B");
	output = findAndReplaceAll(output,"[T]","W");
}

string findAndReplaceAll(string input, const string &searchFor, const string &replaceWith)
{
	string stringOne = input;
	string stringTwo;
	while (1)
	{
		stringTwo = findAndReplaceOnceThrough(stringOne, searchFor,replaceWith);
		stringOne = findAndReplaceOnceThrough(stringTwo, searchFor,replaceWith);
		if (stringOne==stringTwo)
		{
			return stringOne;
		}
	}
}

string findAndReplaceOnceThrough(string input, const string &searchFor, const string &replaceWith)
{
	int start = 0;
	string output;
	for (int i = 0; i < input.size(); i++)
	{
		output = findAndReplace(input, searchFor, replaceWith, start);
		if (output==input)
		{
			return output;
		}
		input = output;
		start += replaceWith.size();
		i = start;
	}
	return input;
}

string findAndReplace(string input, const string &searchFor, const string &replaceWith, int &start)
{

	unsigned int first = input.find(searchFor,start);
	if (string::npos == input.find(searchFor,start))
	{
		return input;
	}
	input.erase(first,searchFor.size());
	input.insert(first,replaceWith);
	start = first;
	return input;
}

bool isLoop(const string &one, const string &two)
{
	if (one.size()!=two.size())
	{
		return 0;
	}
	if (RNAtoDNA(one) == DNAantiCodon(flipString(RNAtoDNA(two))))
	{
		return 1;
	}
	return 0;
}

string getTypes(const string &code, const vector<Codon> &codonSheet)
{
	vector<int> temp;
	return getTypes(code, temp, codonSheet);
}

string getTypes(const string &code, vector<int> &types, const vector<Codon> &codonSheet)
{
	ostringstream temp;
	temp << "The code type(s) is/are: ";
	if (isDNA(code))
	{
		temp << "DNA";
		types.push_back(0);
	}
	if (isRNA(code))
	{
		if (types.size())
		{
			temp << ", ";
		}
		temp << "RNA";
		types.push_back(1);
	}
	if (isNucleotide(code))
	{
		if (types.size())
		{
			temp << ", ";
		}
		temp << "Nucleotide";
		types.push_back(2);
	}
	if (isPeptideLetter(code, codonSheet))
	{
		if (types.size())
		{
			temp << ", ";
		}
		temp << "Peptide";
		types.push_back(3);
	}
	if (isPeptideThreeLetters(code, codonSheet))
	{
		if (types.size())
		{
			temp << ", ";
		}
		temp << "Peptide31";
		types.push_back(4);
	}
	return temp.str();
}

bool printToLog(string input)
{
	string fileLocation = "log.txt";
	ofstream file;
	file.open(fileLocation.c_str(), ios::out | ios::app);
	if (file.fail())
	{
		return 0;
	}
	file << input << endl;
	file.close();
	return 1;
}


string threeLetterPeptideToLetterPeptide(const string &peptides, vector<Codon> codonSheet)
{
	ostringstream sequence;

	if(!isPeptideThreeLetters(peptides, codonSheet))
	{
		printToLog("Error: Does not contain all Peptides");
		return "";
	}
	for (int i = 0; i < peptides.size(); i=i+3)
	{
		string tempString;
		tempString.push_back(peptides.at(i));
		tempString.push_back(peptides.at(i+1));
		tempString.push_back(peptides.at(i+2));

		for (int j = 0; j < codonSheet.size(); j++)
		{
			if (tempString==codonSheet.at(j).peptide3l)
			{
				sequence << codonSheet.at(j).peptide;
				break;
			}
		}
	}
	return sequence.str();
}


string letterPeptideToThreeLetterPeptide(const string &peptides, vector<Codon> codonSheet)
{
	ostringstream sequence;
	string tempString;
	if(!isPeptideLetter(peptides, codonSheet))
	{
		printToLog("Error: Does not contain all Peptides");
		return "";
	}
	for (int i = 0; i < peptides.size(); i++)
	{
		for (int j = 0; j < codonSheet.size(); j++)
		{
			if (peptides.at(i)==codonSheet.at(j).peptide)
			{
				sequence << codonSheet.at(j).peptide3l;
				break;
			}
		}
	}
	return sequence.str();
}

bool isDNA(const string &DNA)
{
	for (int i = 0; i < DNA.size(); i++)
	{
		if (!(DNA.at(i)=='A' || DNA.at(i)=='T' || DNA.at(i) =='G' || DNA.at(i) =='C'))
		{
			return 0;
		}
	}
	return 1;
}

bool isRNA(const string &RNA)
{
	for (int i = 0; i < RNA.size(); i++)
	{
		if (!(RNA.at(i)=='A' || RNA.at(i)=='U' || RNA.at(i) =='G' || RNA.at(i) =='C'))
		{
			return 0;
		}
	}
	return 1;
}

bool isNucleotide(const string &bases)
{
	for (int i = 0; i < bases.size(); i++)
	{
		if (!(bases.at(i)=='A' || bases.at(i)=='T' || bases.at(i) =='G' || bases.at(i) =='C' || bases.at(i) =='U'))
		{
			return 0;
		}
	}
	return 1;
}

bool isPeptideLetter(const string &peptides, const vector<Codon> &codonSheet)
{
	for (int i = 0; i < peptides.size(); i++)
	{
		bool good = 0;
		for (int j = 0; j < codonSheet.size(); j++)
		{
			if (peptides.at(i)==codonSheet.at(j).peptide)
			{
				good = 1;
			}
		}
		if (good==0)
		{
			return 0;
		}
	}
	return 1;
}

bool isPeptideThreeLetters(const string &peptides, const vector<Codon> &codonSheet)
{
	if (!(peptides.size()%3==0))
	{
		return 0;
	}
	for (int i = 0; i < peptides.size(); i=i+3)
	{

		string tempString;
		tempString.push_back(i);
		tempString.push_back(i+1);
		tempString.push_back(i+2);
		bool good = 0;
		for (int j = 0; j < codonSheet.size(); j++)
		{
			if (tempString==codonSheet.at(j).peptide3l)
			{
				good = 1;
			}
		}
		if (good==0)
		{
			return 0;
		}
	}
	return 1;
}

vector<Frame> getMinimumFrames(const string &DNA, const vector<Codon> &codonSheet, const unsigned int &minimumSize)
{
	vector<Frame> forwardFrames = getFrames(DNA,codonSheet);
	vector<Frame> reverseFrames = getFrames(DNAflipAntiString(DNA), codonSheet);
	for (int i = 0; i < forwardFrames.size(); i++)
	{
		forwardFrames.at(i).marker = 0;
	}
	for (int i = 0; i < reverseFrames.size(); i++)
	{
		reverseFrames.at(i).marker = 1;
	}
	vector<Frame> Frames = mergeFrames(forwardFrames, reverseFrames);
	vector<Frame> bigFrames = findMinFrames(Frames,minimumSize);
	return bigFrames;
}
/* TO DO

-build counters
	-amino
	-nucleotides
	-codons
-fix spliceSites stuff
-work on folding algorithm
-add prefix to header files
-find and replace
-replace "int" with "unsigned long int"

-openGL, openCL, openMP and Qt

*/

/*
Peptide getPeptideIUPACCodon(Peptide peptide)
{
	for (int i = 0; i < peptide.bases.size(); i++)
	{

		string bases;
		peptide.bases.at(i).threeLetters.at*

	}

}
*/

void addMasters(vector<Peptide> &peptides)
{
	for (int i = 0; i < peptides.size(); i++)
	{

		for (int j = 0; j < 3; j++)
		{
			vector<char> letters;
			for (int k = 0; k < peptides.at(i).bases.size(); k++)
			{
				letters.push_back(peptides.at(i).bases.at(j).threeLetters.at(k));
			}


			//peptides.at(i).IUPACBases = getMasters(letters);
		}

	}


}


vector<char> getMasters(const vector<char> &characters)
{
	char R = 0;
	char Y = 0;
	char K = 0;
	char M = 0;
	char S = 0;
	char W = 0;
	char B = 0;
	char D = 0;
	char H = 0;
	char V = 0;
	char N = 0;
	vector<char> output;
	for (int i = 0; i < characters.size(); i++)
	{
		if (isSubset('R',characters.at(i)))
		{
			R++;
		}
		if (isSubset('Y',characters.at(i)))
		{
			Y++;
		}
		if (isSubset('K',characters.at(i)))
		{
			K++;
		}
		if (isSubset('M',characters.at(i)))
		{
			M++;
		}
		if (isSubset('S',characters.at(i)))
		{
			S++;
		}
		if (isSubset('W',characters.at(i)))
		{
			W++;
		}
		if (isSubset('B',characters.at(i)))
		{
			B++;
		}
		if (isSubset('D',characters.at(i)))
		{
			D++;
		}
		if (isSubset('H',characters.at(i)))
		{
			H++;
		}
		if (isSubset('V',characters.at(i)))
		{
			V++;
		}
		if (isSubset('H',characters.at(i)))
		{
			V++;
		}
	}
	if (R >= 2)
	{
		output.push_back('R');
	}
	if (Y >= 2)
	{
		output.push_back('Y');
	}
	if (K >= 2)
	{
		output.push_back('K');
	}
	if (M >= 2)
	{
		output.push_back('M');
	}
	if (S >= 2)
	{
		output.push_back('S');
	}
	if (W >= 2)
	{
		output.push_back('W');
	}

	if (B >= 3)
	{
		output.push_back('B');
	}
	if (D >= 3)
	{
		output.push_back('D');
	}
	if (H >= 3)
	{
		output.push_back('H');
	}
	if (V >= 3)
	{
		output.push_back('V');
	}
	if (N >= 4)
	{
		output.push_back('N');
	}
	return output;
}

char fromIUPAC(const char &large)
{
	char small;
	if (large == 'T' || large == 'U' || large == 'A' || large == 'C' || large == 'G')
	{
		small = large;
	}
	else if (large == 'R')
	{
		int random = getRandomNumber(0,1,1);
		if (random==0)
		{
			small = 'A';
		}
		else if (random==1)
		{
			small = 'G';
		}
	}
	else if (large == 'Y')
	{
		int random = getRandomNumber(0,1,1);
		if (random==0)
		{
			small = 'C';
		}
		else if (random==1)
		{
			small = 'T';
		}
	}
	else if (large == 'K')
	{
		int random = getRandomNumber(0,1,1);
		if (random==0)
		{
			small = 'G';
		}
		else if (random==1)
		{
			small = 'T';
		}
	}
	else if (large == 'M')
	{
		int random = getRandomNumber(0,1,1);
		if (random==0)
		{
			small = 'A';
		}
		else if (random==1)
		{
			small = 'C';
		}
	}
	else if (large == 'S')
	{
		int random = getRandomNumber(0,1,1);
		if (random==0)
		{
			small = 'C';
		}
		else if (random==1)
		{
			small = 'G';
		}

	}
	else if (large == 'W')
	{
		int random = getRandomNumber(0,1,1);
		if (random==0)
		{
			small = 'A';
		}
		else if (random==1)
		{
			small = 'T';
		}
	}
	else if (large == 'B')
	{
		int random = getRandomNumber(0,2,1);
		if (random==0)
		{
			small = 'G';
		}
		else if (random==1)
		{
			small = 'T';
		}
		else if (random==2)
		{
			small = 'C';
		}
	}
	else if (large == 'D')
	{
		int random = getRandomNumber(0,2,1);
		if (random==0)
		{
			small = 'G';
		}
		else if (random==1)
		{
			small = 'T';
		}
		else if (random==2)
		{
			small = 'A';
		}
	}
	else if (large == 'H')
	{
		int random = getRandomNumber(0,2,1);
		if (random==0)
		{
			small = 'C';
		}
		else if (random==1)
		{
			small = 'T';
		}
		else if (random==2)
		{
			small = 'A';
		}
	}
	else if (large == 'H')
	{
		int random = getRandomNumber(0,2,1);
		if (random==0)
		{
			small = 'C';
		}
		else if (random==1)
		{
			small = 'G';
		}
		else if (random==2)
		{
			small = 'A';
		}
	}
	else if (large == 'N')
	{
		int random = getRandomNumber(0,3,1);
		if (random==0)
		{
			small = 'C';
		}
		else if (random==1)
		{
			small = 'G';
		}
		else if (random==2)
		{
			small = 'A';
		}
		else if (random==2)
		{
			small = 'T';
		}
	}
	else
	{
		small = ' ';
	}
	return small;
}

string fromIUPAC(const string &DNA)
{
	string output;
	for (int i = 0; i < DNA.size(); i++)
	{
		output.push_back(fromIUPAC(DNA.at(i)));
	}
}





ThreeLettersCount pickBase(vector<ThreeLettersCount> bases)
{
	unsigned int totalNumber = 0;
	for (int i = 0; i < bases.size(); i++)
	{
		totalNumber += bases.at(i).number;
	}
	unsigned int randomNumber = getRandomNumber(0, totalNumber, 1);

	unsigned int countUp = 0;

	for (int i = 0; i < bases.size(); i++)
	{
		if (randomNumber < (bases.at(i).number + countUp))
		{
			return bases.at(i);
		}
		countUp += bases.at(i).number;
	}
}

string getSequenceFromProtein(const Frame &proteinFrame, vector<Peptide> &peptides)
{
	ostringstream protein;
	for (int i = 0; i < proteinFrame.list.size(); i++)
	{
		protein << proteinFrame.list.at(i).peptide;
	}
	string proteinString = protein.str();
	//string proteinString = "MKDLNNTKGNTKSEGSTERGNSGVDRGIVVPNTQIKMRFLNQVRYYSVNNNLKIGKDTNIELSKDTSTSDLLEFEKLVIDNINEENINNNLLSIIKNVDILILAYNRIKSKPGNITPGTTLETLDGINIIYLNKLSNELGTGKFKFKPMRIVNIPKPKGGIRPLSVGNPRDKIVQEVIRIILDTIFDKKISTHSHGFRKNISCQTAI";
	ostringstream sequence;
	for (int i = 0; i < proteinString.size(); i++)
	{
		for (int j = 0; j < peptides.size(); j++)
		{
			if (peptides.at(j).codon.peptide == proteinString.at(i))
			{
				int randomNumber = getRandomNumber(0,peptides.at(j).countBases(),1);
				unsigned int countUp = 0;
				ThreeLettersCount basePicked = pickBase(peptides.at(j).bases);
				sequence << basePicked.threeLetters;
				/*
				for (int k = 0; k < peptides.at(j).bases.size(); k++)
				{
					//if ((peptides.at(j).bases.at(k).number + countUp) > randomNumber)
					{
						sequence << peptides.at(j).bases.at(k).threeLetters;
						break;
					}
					//else
					{
						countUp += peptides.at(j).bases.at(k).number;
					}
				}
				*/

			}
		}
	}
	for (int j = 0; j < peptides.size(); j++)
	{
		if (peptides.at(j).codon.peptide == 'n')
		{
			int randomNumber = getRandomNumber(0,peptides.at(j).countBases(),1);
			//unsigned int countUp = 0

			ThreeLettersCount basePicked = pickBase(peptides.at(j).bases);
			sequence << basePicked.threeLetters;
			/*
			for (int k = 0; k < peptides.at(j).bases.size(); k++)
			{
				if ((peptides.at(j).bases.at(k).number + countUp) > randomNumber)
				{
					sequence << peptides.at(j).bases.at(k).threeLetters;
					break;
				}
				else
				{
					countUp += peptides.at(j).bases.at(k).number;
				}
			}
			*/

		}
	}
	return sequence.str();


}

string getSequenceFromProtein(const string proteinString, vector<Peptide> &peptides)
{
	ostringstream sequence;
	for (int i = 0; i < proteinString.size(); i++)
	{
		for (int j = 0; j < peptides.size(); j++)
		{
			if (peptides.at(j).codon.peptide == proteinString.at(i))
			{
				int randomNumber = getRandomNumber(0,peptides.at(j).countBases(),1);
				unsigned int countUp = 0;
				ThreeLettersCount basePicked = pickBase(peptides.at(j).bases);
				sequence << basePicked.threeLetters;
				/*
				for (int k = 0; k < peptides.at(j).bases.size(); k++)
				{
					//if ((peptides.at(j).bases.at(k).number + countUp) > randomNumber)
					{
						sequence << peptides.at(j).bases.at(k).threeLetters;
						break;
					}
					//else
					{
						countUp += peptides.at(j).bases.at(k).number;
					}
				}
				*/

			}
		}
	}
	for (int j = 0; j < peptides.size(); j++)
	{
		if (peptides.at(j).codon.peptide == 'n')
		{
			int randomNumber = getRandomNumber(0,peptides.at(j).countBases(),1);
			//unsigned int countUp = 0

			ThreeLettersCount basePicked = pickBase(peptides.at(j).bases);
			sequence << basePicked.threeLetters;
			/*
			for (int k = 0; k < peptides.at(j).bases.size(); k++)
			{
				if ((peptides.at(j).bases.at(k).number + countUp) > randomNumber)
				{
					sequence << peptides.at(j).bases.at(k).threeLetters;
					break;
				}
				else
				{
					countUp += peptides.at(j).bases.at(k).number;
				}
			}
			*/

		}
	}
	return sequence.str();


}



int getRandomNumber(const int &minimum, const int &maximum, const int &iteration)
{
	return rand()%maximum+minimum;
}

string getStatistics(const vector<Codon> &codonSheet, const string &DNA, bool &error, vector<Peptide> &peptides)
{
	bool usePeptideData = 1;

	ostringstream dataString;
	dataString << "Statistics:" << endl << "-----------" << endl << endl;
	dataString << "Nucleotides:" << endl;
	vector<LetterCount> letters = countLetters(DNA);
	unsigned int totalLetters = 0;
	for (int i = 0; i < letters.size(); i++)
	{
		totalLetters += letters.at(i).number;
	}
	for (int i = 0; i < letters.size(); i++)
	{
		dataString << letters.at(i).letter << "\t" << letters.at(i).number << "\t(" << (double)letters.at(i).number/totalLetters*100 << "%)" << endl;
	}
	dataString << endl;

	dataString << "Codons:" << endl;
	vector<ThreeLettersCount> threeLetters = countThreeLetters(DNA,error);
	int totalThreeLetters = 0;
	for (int i = 0; i < threeLetters.size(); i++)
	{
		totalThreeLetters += threeLetters.at(i).number;
	}
	for (int i = 0; i < threeLetters.size(); i++)
	{
		dataString << threeLetters.at(i).threeLetters << "\t" << threeLetters.at(i).number << "\t(" << (double) threeLetters.at(i).number/totalThreeLetters*100 << "%)" << endl;
	}
	dataString << endl;











	if (usePeptideData)
	{
			dataString << "Peptides:" << endl;
		//vector<Peptide> peptides;
		//Makes peptide sheet.
		for (int i = 0; i < codonSheet.size(); i++)
		{
			bool isThere = 0;
			for (int j = 0; j < peptides.size(); j++)
			{
				if (peptides.at(j).codon.peptide == codonSheet.at(i).peptide && peptides.at(j).codon.peptide3l == codonSheet.at(i).peptide3l)
					// && peptides.at(j).codon.peptideName == codonSheet.at(i).peptideName
				{
					isThere = 1;

					ThreeLettersCount tempThreeLetters;
					tempThreeLetters.threeLetters = codonSheet.at(i).bases;
					tempThreeLetters.number = 0;

					peptides.at(j).bases.push_back(tempThreeLetters);
					break;
				}
			}
			if (isThere == 0)
			{
				Peptide tempPeptide;
				Codon tempCodon = codonSheet.at(i);
				tempCodon.bases = "";
				if (tempCodon.peptide=='n')
				{
					tempCodon.peptideName = "Stop";
				}
				tempPeptide.codon = tempCodon;

				ThreeLettersCount tempThreeLetters;
				tempThreeLetters.threeLetters = codonSheet.at(i).bases;
				tempThreeLetters.number = 0;
				tempPeptide.bases.push_back(tempThreeLetters);

				peptides.push_back(tempPeptide);
			}

		}

	//WORKS ABOVE



		for (int i = 0; i < threeLetters.size(); i++)
		{
			bool foundIt = 0;
			for (int j = 0; j < peptides.size(); j++)
			{
				for (int k = 0; k < peptides.at(j).bases.size(); k++)
				{
					if (threeLetters.at(i).threeLetters == peptides.at(j).bases.at(k).threeLetters)
					{
						peptides.at(j).bases.at(k).number = threeLetters.at(i).number;
						foundIt = 1;
						break;
					}
				}
				if (foundIt == 1)
				{
					break;
				}
			}
		}

		for (int i = 0; i < peptides.size(); i++)
		{
			dataString << peptides.at(i).codon.peptide << "\t" << peptides.at(i).countBases() << "(" << (double)peptides.at(i).countBases()/totalThreeLetters*100 << "%)\t";
			for (int j = 0; j < peptides.at(i).bases.size(); j++)
			{
				dataString << peptides.at(i).bases.at(j).threeLetters << "(";
				if (peptides.at(i).countBases() != 0)
				{
					dataString << (double)peptides.at(i).bases.at(j).number/peptides.at(i).countBases()*100;
				}
				else
				{
					dataString << "NaN";
				}
				dataString << "%)"<< "\t";
			}
			dataString << endl;
		}
	}



	/*
	//Use Statistics to find
	string protein;
	for (int i = 0; i < protein.size(); i++)
	{
		for (int j = 0; j < peptides.size(); j++)
		if protein.at()

	}
	*/


	/*
	for (int i = 0; i < peptides.size(); i++)
	{
		dataString << peptides.at(i).codon.peptide << "\t";
		for (int j = 0; j < peptides.at(i).bases.size(); j++)
		{
			dataString << peptides.at(i).bases.at(j).threeLetters << "(";
			if (peptides.at(i).countBases() != 0)peptides.at(i).bases.at(j).threeLetters
			{
				dataString << peptides.at(i).bases.at(j).number/peptides.at(i).countBases()*100;
			}
			else
			{
				dataString << "NaN";
			}
			dataString << "%)"<< "\t";
		}
		dataString << endl;
	}
	*/
	/*
	for (int j = 0; j < threeLetters.size(); j++)
		{
			if (peptides.at(i).bases.atthreeLetters == threeLetters.at(j).threeLetters)
			{
				Peptide temp;
				temp.codon.
				peptides.at(i).
			}
		}
	dataString << codonSheet.at().bases
	*/

	return dataString.str();
}

//vector<Codon36Count> countThreeLetters(const string &code)

vector<Amino> countPeptide(const vector<Codon> &codons)
{
	vector <Amino> data;
	for(unsigned int i=0; i < codons.size(); i++)
	{
		bool itemFound = false;
		for (unsigned int j=0; j < data.size(); j++)
		{
			if (codons.at(i).peptide==data.at(j).peptide && codons.at(i).peptide3l==data.at(j).peptide3l && codons.at(i).peptideName==data.at(j).peptideName)
			{
				itemFound = true;
				data.at(j).number++;
				break;
			}
		}
		if (itemFound==false)
		{
			Amino temp;
			temp.number = 1;
			temp.peptide = codons.at(i).peptide;
			temp.peptide3l = codons.at(i).peptide3l;
			temp.peptideName = codons.at(i).peptideName;
			data.push_back(temp);
		}
	}
	return data;
}

vector<ThreeLettersCount> countCodons(const vector<Codon> &codons)
{
	vector <ThreeLettersCount> data;
	for(unsigned int i=0; i < codons.size(); i++)
	{
		bool itemFound = false;
		for (unsigned int j=0; j < data.size(); j++)
		{
			if (codons.at(i).bases==data.at(j).threeLetters)
			{
				itemFound = true;
				data.at(j).number++;
				break;
			}
		}
		if (itemFound==false)
		{
			ThreeLettersCount temp;
			temp.number = 1;
			temp.threeLetters = codons.at(i).bases;
			data.push_back(temp);
		}
	}
	return data;
}

//Probably will never be used.
vector<ThreeLettersCount> countThreeLetters(const string &code, bool &error)
{
	vector <ThreeLettersCount> data;
	for (int i = 0; i < code.size() - 2; i = i + 3)
	{
		string tmpLetters;
		tmpLetters.clear();
		if (i + 2 >= code.size() || error == 1)
		{
			error = 1;
			break;
		}
		tmpLetters.push_back(code.at(i));
		tmpLetters.push_back(code.at(i + 1));
		tmpLetters.push_back(code.at(i + 2));
		//Codon tmpLetters = findCodon(tmpBases, codonSheet);

		//if (tmpCodon.peptideName != "")
		//{
			//list.push_back(tmpCodon);
		//
		bool itemFound = false;
		for (unsigned int j=0; j < data.size(); j++)
		{
			if (tmpLetters==data.at(j).threeLetters)
			{
				itemFound = true;
				data.at(j).number++;
				break;
			}
		}
		if (itemFound==false)
		{
			ThreeLettersCount temp;
			temp.number = 1;
			temp.threeLetters = tmpLetters;
			data.push_back(temp);
		}
	}
	return data;
}


vector<LetterCount> countLetters(const string &code)
{
	vector <LetterCount> data;
	for(unsigned int i=0; i < code.size(); i++)
	{
		bool itemFound = false;
		for (unsigned int j=0; j < data.size(); j++)
		{
			if (code.at(i)==data.at(j).letter)
			{
				itemFound = true;
				data.at(j).number++;
				break;
			}
		}
		if (itemFound==false)
		{
			LetterCount temp;
			temp.number = 1;
			temp.letter = code.at(i);
			data.push_back(temp);
		}
	}
	return data;
}

double measureDistance(Point &point1, Point &point2)
{
	double distance;
	distance = sqrt(pow(point1.x-point2.x,2)+pow(point1.y-point2.y,2)+pow(point1.z-point2.z,2));
	return distance;
}

Point translate(const Point &input, const double &x, const double &y, const double &z)
{
	Point output;
	output.x=input.x+x;
	output.y=input.y+y;
	output.z=input.z+z;
	return output;
}

Point stretch(const Point &input, const double &x, const double &y, const double &z)
{
	Point output;
	output.x=input.x*x;
	output.y=input.y*y;
	output.z=input.z*z;
	return output;
}

Point rotate_x(const Point &input, const double &angle)
{
	Point output;
	output.x=input.x;
	output.y=input.y*cos(angle)-input.z*sin(angle);
	output.z=input.z*cos(angle)+input.y*sin(angle);
	return output;
}

Point rotate_y(const Point &input, const double &angle)
{
	Point output;
	output.x=input.x*cos(angle)+input.z*sin(angle);
	output.y=input.y;
	output.z=input.z*cos(angle)-input.x*sin(angle);
	return output;
}

Point rotate_z(const Point &input, const double &angle)
{
	Point output;
	output.x=input.x*cos(angle)-input.y*sin(angle);
	output.y=input.y*cos(angle)+input.x*sin(angle);
	output.z=input.z;
	return output;
}

double measureAngle(Point &point1, Point &point2, Point &point3)
{
	double distance_A = measureDistance(point1, point2);
	double x_A = point2.x-point1.x;
	double y_A = point2.y-point1.y;
	double z_A = point2.z-point1.z;

	cout << endl << x_A << " " << y_A << " " << z_A << " " << distance_A;

	double distance_B = measureDistance(point3, point2);
	double x_B = point2.x-point3.x;
	double y_B = point2.y-point3.y;
	double z_B = point2.z-point3.z;

	cout << endl << x_B << " " << y_B << " " << z_B << " " << distance_B;

	double dotProduct;
	dotProduct = x_A*x_B + y_A*y_B + z_A*z_B;

	cout << endl << dotProduct << endl;

	double cosineOfAngle = dotProduct/(distance_A*distance_B);
	cout << endl << cosineOfAngle << endl;

	if (cosineOfAngle < -1)
	{
		cosineOfAngle = -1;
	}
	if (cosineOfAngle > 1)
	{
		cosineOfAngle = 1;
	}
	double angle = acos(cosineOfAngle); //In radians
	//angle = angle*180/3.14159265359; //Converts to Degrees
	return angle;
}


/*
//Not sure if this works!
string findSpliceString(const string &code, const SpliceSite &splicePoints, int startHere, int &startNum, int &firstNum, int &lastNum, int &stopNum)
{
	string startFirst = splicePoints.start + splicePoints.first;
	string lastStop = splicePoints.last + splicePoints.stop;
	startNum = code.find(startFirst,startHere);
	firstNum = code.find(startFirst,startNum+splicePoints.start.size());

	firstNum+splicePoints.first.size();
	middleWorks = true;
	for (int i=0; i < splicePoints.middle.size(); i++)
	{

	}
	//insert breaks

	last = code.find(,startFirst.size()+startHere);


}
*/
/*
string toIUPAC(const string &code)
{
	findAndReplace(code, "{A}", "B");
	findAndReplace(code, "{C}", "D");
	findAndReplace(code, "{G}", "H");
	findAndReplace(code, "{T}", "V");
	findAndReplace(code, "{T}", "V");
	findAndReplace(code, "{T}", "V");
	findAndReplace(code, "{T}", "V");
	findAndReplace(code, "{T}", "V");
}
*/

/*
string findAndReplace(const string &main, const string &find, const string &replace)
{
	int startSearch = 0;
	vector<int> locations;
	for (startSearch = 0; startSearch < main.size(); startSearch++)
	{
		int location = main.find(find, startSearch);
		main.replace(location, find.size());//location, find.size());
		//main.insert(location, replace);
	}

}
*/

bool isSubset(const char &large, const char &small)
{
	bool subset = false;
	if (large == small)
	{
		subset = true;
	}
	else if ((small == 'T' || small == 'U') && (large == 'T' || large == 'U'))
	{
		subset = true;
	}
	else if (small == 'C' && (large == 'S' || large == 'M' || large == 'Y' || large == 'B' || large == 'H' || large == 'V' || large == 'N'))
	{
		subset = true;
	}
	else if (small == 'A' && (large == 'W' || large == 'M' || large == 'R' || large == 'D' || large == 'H' || large == 'V' || large == 'N'))
	{
		subset = true;
	}
	else if (small == 'G' && (large == 'S' || large == 'K' || large == 'R' || large == 'B' || large == 'D' || large == 'V' || large == 'N'))
	{
		subset = true;
	}
	else if ((small == 'T' || small == 'U') && (large == 'W' || large == 'K' || large == 'Y' || large == 'B' || large == 'D' || large == 'H' || large == 'N'))
	{
		subset = true;
	}
	else if (small == 'W' && (large == 'D' || large == 'H' || large == 'N'))
	{
		subset = true;
	}
	else if (small == 'S' && (large == 'B' || large == 'V' || large == 'N'))
	{
		subset = true;
	}
	else if (small == 'M' && (large == 'H' || large == 'V' || large == 'N'))
	{
		subset = true;
	}
	else if (small == 'K' && (large == 'B' || large == 'D' || large == 'N'))
	{
		subset = true;
	}
	else if (small == 'R' && (large == 'D' || large == 'V' || large == 'N'))
	{
		subset = true;
	}
	else if (small == 'Y' && (large == 'B' || large == 'H' || large == 'N'))
	{
		subset = true;
	}
	else if (small == 'B' && large == 'N')
	{
		subset = true;
	}
	else if (small == 'D' && large == 'N')
	{
		subset = true;
	}
	else if (small == 'H' && large == 'N')
	{
		subset = true;
	}
	else if (small == 'V' && large == 'N')
	{
		subset = true;
	}
	return subset;
}

vector<Frame> mergeFrames(const vector<Frame> &Frames1, const vector<Frame> &Frames2)
{
	vector<Frame> mergedFrames;
	for (int i = 0; i < Frames1.size(); i++)
	{
		mergedFrames.push_back(Frames1.at(i));
	}
	for (int i = 0; i < Frames2.size(); i++)
	{
		mergedFrames.push_back(Frames2.at(i));
	}
	return mergedFrames;
}

vector<Frame> findMinFrames(const vector<Frame> &Frames, const int &minsize)
{
	vector<Frame> bigFrames;
	for (int i = 0; i < Frames.size(); i++)
	{
		if (Frames.at(i).list.size() > minsize)
		{
			bigFrames.push_back(Frames.at(i));
		}
	}
	return bigFrames;
}

vector<Paladrome> findPaladromes(const string &code, const int &minsize)
{
	vector<Paladrome> Paladromes;
	for (int i = 0; i < code.size(); i++)
	{
		string tmpString;
		for (int j = i; j < code.size(); j++)
		{
			//if (j >= code.size())
			//{
				//break;
			//}
			tmpString.push_back(code.at(j));
			//cout << endl << tmpString << endl;
			if (isPaladrome(tmpString)&&(tmpString.size()>=minsize))
			{
				Paladrome tmpPaladrome;
				tmpPaladrome.sequence = tmpString;
				tmpPaladrome.start = i;
				tmpPaladrome.stop = j;
				Paladromes.push_back(tmpPaladrome);
			}
		}
		tmpString.clear();
	}
	return Paladromes;
}

string DNAflipAntiString(string DNA)
{
	string flipedAntiString = flipString(DNAantiCodon(RNAtoDNA(DNA)));
	return flipedAntiString;
}

bool isPaladrome(const string Code)
{
	//Perhaps add IUPAC. Already did>
	if (flipString(RNAtoDNA(Code)) == DNAantiCodon(RNAtoDNA(Code)))
	{
		return true;
	}
	return false;
}

bool printCodonSheetToFile(const vector<Codon> &codonSheet, const string &filePath)
{
	ofstream fout;
	bool error = 0;
	fout.open(filePath.c_str());
	if (fout.fail())
	{
		printToLog("Error: Failed to print the codon sheet.");
		//cout << "Failed to print the codon sheet." << endl;
		error = 1;
	}
	fout << "DNA Codon	Peptide Letter	Peptide Name	Peptide 3L	Special	Start	Stop	\"Acidic = 1, Basic = 2, Polar = 3, Non - Polar = 4\"	Bond Length	Bond Angle" << endl;
	for (int i = 0; i < codonSheet.size(); i++)
	{
		fout << codonSheet.at(i).bases << "\t" << codonSheet.at(i).peptide << "\t" << codonSheet.at(i).peptideName << "\t" << codonSheet.at(i).peptide3l << "\t" << codonSheet.at(i).special << "\t" << codonSheet.at(i).start << "\t" << codonSheet.at(i).stop << "\t" << codonSheet.at(i).type << "\t" << codonSheet.at(i).bondLength << "\t" << codonSheet.at(i).bondAngle << "\t" << endl;
	}
	//cout << "The size of the Codon Sheet is: " << codonSheet.size() << endl;
	fout.close();
	return error;
}

vector<Codon> getCodonSequence(int startCodon, int &lastCodon, const vector<Codon> &codonSheet, const string &DNA, string &codingDNA, bool &error, bool &doesBreak)
{
	vector<Codon> list;
	if (startCodon == -1)
	{
		error = 1;
	}
	for (int i = startCodon; i < DNA.size() - 2; i = i + 3)
	{
		string tmpBases;
		tmpBases.clear();
		if (i + 2 >= DNA.size() || error == 1)
		{
			error = 1;
			break;
		}
		tmpBases.push_back(DNA.at(i));
		tmpBases.push_back(DNA.at(i + 1));
		tmpBases.push_back(DNA.at(i + 2));
		Codon tmpCodon = findCodon(tmpBases, codonSheet);
		if (tmpCodon.stop == 1)
		{
			codingDNA.append(tmpCodon.bases);
			doesBreak = 1;
			break;
		}
		//if (tmpCodon.peptideName != "")
		//{
			list.push_back(tmpCodon);
		//}
		codingDNA.append(tmpBases);
		lastCodon = i + 5;
	}
	/*
	for (int i = startCodon; i <= lastCodon; i++)
	{
		codingDNA.push_back(DNA.at(i));
	}
	*/
	return list;
}

//Test this
vector<Codon> getAminoSequence(const vector<Codon> &codonSheet, const string &amino, bool &error)
{
	vector<Codon> list;
	//findCodon(tmpBases, codonSheet);
	for (int i = 0; i < amino.size(); i++)
	{
		Codon tmpCodon;
		for (int j = 0; j < codonSheet.size(); j++)
		{
			if (codonSheet.at(j).peptide == amino.at(j))
			{
				tmpCodon = codonSheet.at(j);
				break;
			}
			tmpCodon.peptide = 'a';
			tmpCodon.peptide3l = "aaa";
			tmpCodon.start = 0;
			tmpCodon.stop = 0;
			tmpCodon.bases = "aaa";
		}
		list.push_back(tmpCodon);

	}
	return list;
}


bool printCodeToFile(string Code, string filePath)
{
	bool error = 0;
	ofstream fout;
	fout.open(filePath.c_str());
	if (fout.fail())
	{
		ostringstream temp;
		temp << "Error: Output stream to \"" << filePath << "\" failed to open";
		printToLog(temp.str());
		//cout << "Output stream to \"" << filePath << "\" failed to open";
		error = 1;
	}
	fout << Code << endl;
	fout.close();
	return error;
}

void printCodeToConsole(string Code)
{
	cout << Code << endl;
}

void printFramesToConsole(const vector<Frame> &Frames)
{
	for (int j = 0; j < Frames.size(); j++)
	{
		for (int i = 0; i < Frames.at(j).list.size(); i++)
		{
			cout << Frames.at(j).list.at(i).peptide;
		}
		cout << endl << endl;

		for (int i = 0; i < Frames.at(j).codingDNA.size(); i++)
		{
			cout << Frames.at(j).codingDNA.at(i);
		}
		cout << " (" << Frames.at(j).startCodon << ", " << Frames.at(j).lastCodon << ")" << endl << endl << endl;
	}
}

bool printFramesToFile(const vector<Frame> &Frames, const string &filePath)
{
	ofstream fout;
	bool error = 0;
	fout.open(filePath.c_str());
	if (fout.fail())
	{
		ostringstream temp;
		temp << "Error: Output stream to \"" << filePath << "\" failed to open";
		printToLog(temp.str());
		//cout << "Output stream to \"" << filePath << "\" failed to open";
		error = 1;
	}
	for (int j = 0; j < Frames.size(); j++)
	{
		for (int i = 0; i < Frames.at(j).list.size(); i++)
		{
			fout << Frames.at(j).list.at(i).peptide;
		}
		fout << endl << endl;

		for (int i = 0; i < Frames.at(j).codingDNA.size(); i++)
		{
			fout << Frames.at(j).codingDNA.at(i);
		}

		fout << " (" << Frames.at(j).startCodon << ", " << Frames.at(j).lastCodon << ")" << endl << endl << endl;
	}

	fout.close();
	return error;
}

bool printMachineFramesToFile(const vector<Frame> &Frames, const string &filePath)
{
	ofstream fout;
	bool error = 0;
	fout.open(filePath.c_str());
	if (fout.fail())
	{
		ostringstream temp;
		temp << "Error: Output stream to \"" << filePath << "\" failed to open";
		printToLog(temp.str());
		//cout << "Output stream to \"" << filePath << "\" failed to open";
		error = 1;
	}
	fout << "Peptides\tDNA\tStart\tStop" << endl;
	for (int j = 0; j < Frames.size(); j++)
	{

		for (int i = 0; i < Frames.at(j).list.size(); i++)
		{
			fout << Frames.at(j).list.at(i).peptide;
		}
		fout << "\t";

		for (int i = 0; i < Frames.at(j).codingDNA.size(); i++)
		{
			fout << Frames.at(j).codingDNA.at(i);
		}

		fout << "\t" << Frames.at(j).startCodon << "\t" << Frames.at(j).lastCodon << endl;
	}

	fout.close();
	return error;
}

vector<Frame> getFrames(const string &DNA, const vector<Codon> &codonSheet)
{
	vector<Frame> Frames;
	for (int j = 0; j < DNA.size();)
	{
		int startCodon = findStartCodon(DNA, codonSheet, j);
		if (startCodon == -1)
		{
			break;
		}

		string codingDNA;
		bool doesBreak = 0;
		bool error = 0;
		int lastCodon;
		vector<Codon> list = getCodonSequence(startCodon, lastCodon, codonSheet, DNA, codingDNA, error, doesBreak);
		Frame tmpFrame;
		tmpFrame.codingDNA = codingDNA;
		tmpFrame.startCodon = startCodon;
		tmpFrame.lastCodon = lastCodon;
		tmpFrame.error = error;
		tmpFrame.doesBreak = doesBreak;
		tmpFrame.list = list;
		Frames.push_back(tmpFrame);
		j = startCodon + 1;
	}
	return Frames;
}

void printCodonSheetToConsole(const vector<Codon> &codonSheet)
{
	cout << "DNA Codon	Peptide Letter	Peptide Name	Peptide 3L	Special	Start	Stop	\"Acidic = 1, Basic = 2, Polar = 3, Non - Polar = 4\"	Bond Length	Bond Angle" << endl;
	for (int i = 0; i < codonSheet.size(); i++)
	{
		cout << codonSheet.at(i).bases << "\t" << codonSheet.at(i).peptide << "\t" << codonSheet.at(i).peptideName << "\t" << codonSheet.at(i).peptide3l << "\t" << codonSheet.at(i).special << "\t" << codonSheet.at(i).start << "\t" << codonSheet.at(i).stop << "\t" << codonSheet.at(i).type << "\t" << codonSheet.at(i).bondLength << "\t" << codonSheet.at(i).bondAngle << "\t" << endl;
	}
	cout << "The size of the Codon Sheet is: " << codonSheet.size() << endl;
}

Codon findCodon(const string &DNA, const vector<Codon> &codonSheet)
{
	for (int i = 0; i < codonSheet.size(); i++)
	{
		if (DNA == codonSheet.at(i).bases)
		{
			return codonSheet.at(i);
		}
	}
	Codon tmp;
	tmp.peptide3l = "aaa";
	tmp.peptide = 'a';
	tmp.start = 0;
	tmp.stop = 0;
	tmp.bases = "aaa";
	return tmp;
}

int findStartCodon(const string &DNA, const vector<Codon> &codonSheet, const int &start)
{
	if (DNA.size() < 2)
	{
		return -1;
	}
	for (int i = start; i < DNA.size() - 2; i++)
	{
		string tmp;
		tmp.push_back(DNA.at(i));
		tmp.push_back(DNA.at(i + 1));
		tmp.push_back(DNA.at(i + 2));
		for (int j = 0; j < codonSheet.size(); j++)
		{
			if (codonSheet.at(j).bases == tmp && codonSheet.at(j).start)
			{
				return i;
			}
		}
	}
	return -1;
}

string RNAtoDNA(string RNA)
{
	for (int i = 0; i < RNA.size(); i++)
	{
		if (RNA.at(i) == 'U')
		{
			RNA.at(i) = 'T';
		}
	}
	return RNA;
}

string DNAtoRNA(string DNA)
{
	for (int i = 0; i < DNA.size(); i++)
	{
		if (DNA.at(i) == 'T')
		{
			DNA.at(i) = 'U';
		}
	}
	return DNA;
}

string DNAantiCodon(string DNA)
{
	bool IUPAC = 0;
	bool basic = 1;
	for (int i = 0; i < DNA.size(); i++)
	{
		if (DNA.at(i) == 'A' && IUPAC == false)
		{
			DNA.at(i) = 'T';
		}
		else if ((DNA.at(i) == 'T' || DNA.at(i) == 'U') && IUPAC == false)
		{
			DNA.at(i) = 'A';
		}
		else if (DNA.at(i) == 'C' && IUPAC == false)
		{
			DNA.at(i) = 'G';
		}
		else if (DNA.at(i) == 'G' && IUPAC == false)
		{
			DNA.at(i) = 'C';
		}
		else if (DNA.at(i) == 'W')
		{
			DNA.at(i) = 'S';
			basic = 0;
		}
		else if (DNA.at(i) == 'S')
		{
			DNA.at(i) = 'W';
			basic = 0;
		}
		else if (DNA.at(i) == 'M')
		{
			DNA.at(i) = 'K';
			basic = 0;
		}
		else if (DNA.at(i) == 'K')
		{
			DNA.at(i) = 'M';
			basic = 0;
		}
		else if (DNA.at(i) == 'R')
		{
			DNA.at(i) = 'Y';
			basic = 0;
		}
		else if (DNA.at(i) == 'Y')
		{
			DNA.at(i) = 'R';
			basic = 0;
		}
		else if (DNA.at(i) == 'B')
		{
			DNA.at(i) = 'A';
			basic = 0;
		}
		else if (DNA.at(i) == 'D')
		{
			DNA.at(i) = 'C';
			basic = 0;
		}
		else if (DNA.at(i) == 'H')
		{
			DNA.at(i) = 'G';
			basic = 0;
		}
		else if (DNA.at(i) == 'V')
		{
			DNA.at(i) = 'T';
			basic = 0;
		}
		else if (DNA.at(i) == 'A')
		{
			DNA.at(i) = 'B';
			basic = 0;
		}
		else if (DNA.at(i) == 'C')
		{
			DNA.at(i) = 'D';
			basic = 0;
		}
		else if (DNA.at(i) == 'G')
		{
			DNA.at(i) = 'H';
			basic = 0;
		}
		else if (DNA.at(i) == 'T')
		{
			DNA.at(i) = 'V';
			basic = 0;
		}
	}
	return DNA;
}

string RNAantiCodon(string RNA)
{
	bool IUPAC = 0;
	bool basic = 1;
	for (int i = 0; i < RNA.size(); i++)
	{
		if (RNA.at(i) == 'A' && IUPAC == false)
		{
			RNA.at(i) = 'U';
		}
		else if ((RNA.at(i) == 'T' || RNA.at(i) == 'U') && IUPAC == false)
		{
			RNA.at(i) = 'A';
		}
		else if (RNA.at(i) == 'C' && IUPAC == false)
		{
			RNA.at(i) = 'G';
		}
		else if (RNA.at(i) == 'G' && IUPAC == false)
		{
			RNA.at(i) = 'C';
		}
		else if (RNA.at(i) == 'W')
		{
			RNA.at(i) = 'S';
			basic = 0;
		}
		else if (RNA.at(i) == 'S')
		{
			RNA.at(i) = 'W';
			basic = 0;
		}
		else if (RNA.at(i) == 'M')
		{
			RNA.at(i) = 'K';
			basic = 0;
		}
		else if (RNA.at(i) == 'K')
		{
			RNA.at(i) = 'M';
			basic = 0;
		}
		else if (RNA.at(i) == 'R')
		{
			RNA.at(i) = 'Y';
			basic = 0;
		}
		else if (RNA.at(i) == 'Y')
		{
			RNA.at(i) = 'R';
			basic = 0;
		}
		else if (RNA.at(i) == 'B')
		{
			RNA.at(i) = 'A';
			basic = 0;
		}
		else if (RNA.at(i) == 'D')
		{
			RNA.at(i) = 'C';
			basic = 0;
		}
		else if (RNA.at(i) == 'H')
		{
			RNA.at(i) = 'G';
			basic = 0;
		}
		else if (RNA.at(i) == 'V')
		{
			RNA.at(i) = 'T';
			basic = 0;
		}
		else if (RNA.at(i) == 'A')
		{
			RNA.at(i) = 'B';
			basic = 0;
		}
		else if (RNA.at(i) == 'C')
		{
			RNA.at(i) = 'D';
			basic = 0;
		}
		else if (RNA.at(i) == 'G')
		{
			RNA.at(i) = 'H';
			basic = 0;
		}
		else if (RNA.at(i) == 'T')
		{
			RNA.at(i) = 'V';
			basic = 0;
		}
	}
	return RNA;
}


string flipString(string input)
{
	for (int i = 0; i < input.size() / 2; i++)
	{
		char tmpChar;
		tmpChar = input.at(i);
		input.at(i) = input.at(input.size() - 1 - i);
		input.at(input.size() - 1 - i) = tmpChar;
	}
	return input;
}

vector<Codon> getCodonSheet(string filePath)
{
	ifstream fin;
	vector<Codon> codonSheet;
	fin.open(filePath.c_str());
	if (fin.fail())
	{
		printToLog("Error: Failed to get the codon sheet.");
		//cout << "Failed to get the codon sheet." << endl;
		return codonSheet;
	}
	string codonLabel;
	getline(fin, codonLabel);
	int i = 0;
	while (1)
	{
		string bases, peptide3l, peptideName;
		char peptide, type;
		double bondLength, bondAngle;
		bool special, start, stop;
		fin >> bases >> peptide >> peptideName >> peptide3l >> special >> start >> stop >> type >> bondLength >> bondAngle;
		Codon tmp;
		tmp.bases = bases;
		tmp.peptide = peptide;
		tmp.peptideName = peptideName;
		tmp.peptide3l = peptide3l;
		tmp.special = special;
		tmp.start = start;
		tmp.stop = stop;
		tmp.type = type;
		tmp.bondLength = bondLength;
		tmp.bondAngle = bondAngle;
		//if (peptideName != "")
		//{
			codonSheet.push_back(tmp);
		//}
		if (fin.eof())
		{
			break;
		}
	}
	fin.close();
	codonSheet.pop_back();
	return codonSheet;
}

string getDNA(string filePath)
{
	ifstream fin;
	fin.open(filePath.c_str());
	if (fin.fail())
	{
		printToLog("Error: Failed to get the DNA.");
		//cout << "Failed to get the DNA." << endl;
		return "0";
	}
	string DNA;
	while (!fin.eof())
	{
		string tmp;
		fin >> tmp;
		//codonSheet.push_back('\n');
		DNA.append(tmp);
	}
	fin.close();
	return DNA;
}
vector<string> parseStrings(string codonString)
{
	string tmpString;
	vector<string> parsedStrings;
	for (int i = 0; i < codonString.size(); i++)
	{
		if (codonString.at(i) == '\n' || codonString.at(i) == ' ' || codonString.at(i) == '\t')
		{
			//string checkWhite = parsedStrings.at(parsedStrings.size() - 1);
			//if (checkWhite != '\n' && checkWhite != ' ' && checkWhite != '\t')
			if (tmpString != "")
			{
				parsedStrings.push_back(tmpString);
			}
			tmpString.clear();
		}
		else if (i == codonString.size() - 1)
		{
			tmpString.push_back(codonString.at(i));
			parsedStrings.push_back(tmpString);
		}
		else
		{
			tmpString.push_back(codonString.at(i));
		}
		return parsedStrings;
	}
}



/*
vector<Codon> parseCodons(string codonString)
{
	vector<Codon> parsedCodons;
	Codon tmpCodon;
	vector<string> parsedStrings = parseStrings(codonString);

	//START HERE WHEN ADDING
	//for (int i = 0; parsedStrings)





}
*/
bool write()
{
	return 0;
}
/*
Codon getAnti(string inputCodon, vector<Codon> codonSheet)
{
	for (int i = 0; i < codonSheet.size(); i++)
	{
		if (codonSheet.at(i).bases == inputCodon)
		{
			return codonSheet.at(i);
		}
	}
	Codon empty;
	return empty;
}

Codon getAnti(string inputCodon, vector<Codon> codonSheet)
{
	for (int i = 0; i < codonSheet.size(); i++)
	{
		if (codonSheet.at(i).bases == inputCodon)
		{
			return codonSheet.at(i);
		}
	}
	Codon empty;
	return empty;
}
*/

/*
for (int i = 0; i < Frames.size(); i++)
{
if (Frames.at(i).list.size() > 20)
{
for (int j = 0; j < Frames.at(i).list.size(); j++)
{
cout << Frames.at(i).list.at(j).peptide;
}
cout << endl << endl;
}

}
*/
