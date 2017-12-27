#include "stdafx.h"
#include <iostream>

#define BLOSUM_SIZE 26
#define MAX_LINE_LENGTH 1024
#define MAX_ARRAY_SIZE 50
#define NEW_LINE_EVERY_CHARS 50
#define FASTA_LINE_MAX_LENGTH 1024

char aminoArray[] = { 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','O' };

void GetBlosum(char* filename, int* resptr);
void ClearString(char* stringPtr, int length);
void ChangeDelim(char* stringPtr, char delim, int length);
int FindInitials(char* stringPtr, char delim, int length, int* resArray, int maxSize);
void PrintBlosum(int* blosum);

void PrintMatrix(int* matrix, int width, int height);
void NeedlemanWunsch(char* seq1, char* seq2, int penalty, int* scoreMat);
void FindStartPoint(int* matrix, int width, int height, int* startWidth, int* startHeight);
void CalMatrices(int matHeight, int matWidth, int penalty, int* gradeMatrix, int* dirMatrix, int* scoreMat, char* seq1, char* seq2);
void FindPath(int * dirMatrix, int startWidth, int startHeight, int matWidth, int matHeight, char* pathPtr);
void GetAlignment(char* path, const char* seq1, const char* seq2, char* resChar1, char* resChar2);
void ReverseString(char* theString);
void PrintResultResp(const char* seq1, const char* seq2, const char* resSeq1, const char* resSeq2);
void PrintResultComp(const char* seq1, const char* seq2, const char* resSeq1, const char* resSeq2);

char* ReadFasta(char* filename);

int main(int argc, char* argv[])
{
	int blosum[BLOSUM_SIZE * BLOSUM_SIZE] = { 4 ,-2 ,0 ,-2 ,-1 ,-2 ,0 ,-2 ,-1 ,0 ,-1 ,-1 ,-1 ,-2 ,-4 ,-1 ,-1 ,-1 ,1 ,0 ,0 ,0 ,-3 ,0 ,-2 ,-1 ,-2 ,4 ,-3 ,4 ,1 ,-3 ,-1 ,0 ,-3 ,0 ,0 ,-4 ,-3 ,3 ,-4 ,-2 ,0 ,-1 ,0 ,-1 ,0 ,-3 ,-4 ,-1 ,-3 ,1 ,0 ,-3 ,9 ,-3 ,-4 ,-2 ,-3 ,-3 ,-1 ,0 ,-3 ,-1 ,-1 ,-2 ,-4 ,-3 ,-3 ,-3 ,-1 ,-1 ,0 ,-1 ,-2 ,-2 ,-2 ,-3 ,-2 ,4 ,-3 ,6 ,2 ,-3 ,-1 ,-1 ,-3 ,0 ,-1 ,-4 ,-3 ,0 ,-4 ,-1 ,0 ,-2 ,0 ,-1 ,0 ,-3 ,-4 ,-1 ,-3 ,1 ,-1 ,1 ,-4 ,2 ,5 ,-3 ,-2 ,0 ,-3 ,0 ,1 ,-3 ,-2 ,-2 ,-4 ,-1 ,2 ,0 ,0 ,-1 ,0 ,-2 ,-3 ,-1 ,-2 ,4 ,-2 ,-3 ,-2 ,-3 ,-3 ,6 ,-3 ,-1 ,0 ,0 ,-3 ,0 ,0 ,-3 ,-4 ,-4 ,-3 ,-3 ,-2 ,-2 ,0 ,-1 ,1 ,-1 ,3 ,-3 ,0 ,-1 ,-3 ,-1 ,-2 ,-3 ,6 ,-2 ,-4 ,0 ,-2 ,-4 ,-3 ,0 ,-4 ,-2 ,-2 ,-2 ,0 ,-2 ,0 ,-3 ,-2 ,-1 ,-3 ,-2 ,-2 ,0 ,-3 ,-1 ,0 ,-1 ,-2 ,8 ,-3 ,0 ,-1 ,-3 ,-2 ,-1 ,-4 ,-2 ,0 ,0 ,-1 ,-2 ,0 ,-3 ,-2 ,-1 ,2 ,0 ,-1 ,-3 ,-1 ,-3 ,-3 ,0 ,-4 ,-3 ,4 ,0 ,-3 ,2 ,1 ,-1 ,-4 ,-3 ,-3 ,-3 ,-2 ,-1 ,0 ,3 ,-3 ,-1 ,-1 ,-3 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,-1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,-1 ,0 ,-3 ,-1 ,1 ,-3 ,-2 ,-1 ,-3 ,0 ,5 ,-2 ,-1 ,0 ,-4 ,-1 ,1 ,2 ,0 ,-1 ,0 ,-2 ,-3 ,-1 ,-2 ,1 ,-1 ,-4 ,-1 ,-4 ,-3 ,0 ,-4 ,-3 ,2 ,0 ,-2 ,4 ,2 ,-4 ,-4 ,-3 ,-2 ,-2 ,-2 ,-1 ,0 ,1 ,-2 ,-1 ,-1 ,-3 ,-1 ,-3 ,-1 ,-3 ,-2 ,0 ,-3 ,-2 ,1 ,0 ,-1 ,2 ,5 ,-1 ,-4 ,-2 ,0 ,-1 ,-1 ,-1 ,0 ,1 ,-1 ,-1 ,-1 ,-1 ,-2 ,3 ,-3 ,1 ,0 ,-3 ,0 ,1 ,-3 ,0 ,0 ,-3 ,-2 ,-1 ,-4 ,-2 ,0 ,0 ,1 ,0 ,0 ,-3 ,-4 ,-1 ,-2 ,0 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,0 ,-4 ,-4 ,-4 ,-4 ,1 ,-4 ,-4 ,-4 ,-4 ,-4 ,0 ,-4 ,-4 ,-4 ,-4 ,-4 ,-1 ,-2 ,-3 ,-1 ,-1 ,-4 ,-2 ,-2 ,-3 ,0 ,-1 ,-3 ,-2 ,1 ,-4 ,7 ,-1 ,-2 ,-1 ,-1 ,0 ,-2 ,-4 ,-2 ,-3 ,-1 ,-1 ,0 ,-3 ,0 ,2 ,-3 ,-2 ,0 ,-3 ,0 ,1 ,-2 ,0 ,0 ,-4 ,-1 ,5 ,1 ,0 ,-1 ,0 ,-2 ,-2 ,-1 ,-1 ,3 ,-1 ,-1 ,-3 ,-2 ,0 ,-3 ,-2 ,0 ,-3 ,0 ,2 ,-2 ,-1 ,0 ,-4 ,-2 ,1 ,5 ,-1 ,-1 ,0 ,-3 ,-3 ,-1 ,-2 ,0 ,1 ,0 ,-1 ,0 ,0 ,-2 ,0 ,-1 ,-2 ,0 ,0 ,-2 ,-1 ,0 ,-4 ,-1 ,0 ,-1 ,4 ,1 ,0 ,-2 ,-3 ,0 ,-2 ,0 ,0 ,-1 ,-1 ,-1 ,-1 ,-2 ,-2 ,-2 ,-1 ,0 ,-1 ,-1 ,-1 ,-3 ,-4 ,-1 ,-1 ,-1 ,1 ,5 ,0 ,0 ,-2 ,0 ,-2 ,-1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,-3 ,-1 ,-3 ,-2 ,-1 ,-3 ,-3 ,3 ,0 ,-2 ,1 ,1 ,-3 ,-4 ,-2 ,-2 ,-3 ,-2 ,0 ,0 ,4 ,-3 ,-1 ,-1 ,-2 ,-3 ,-4 ,-2 ,-4 ,-3 ,1 ,-2 ,-2 ,-3 ,0 ,-3 ,-2 ,-1 ,-1 ,-4 ,-4 ,-2 ,-3 ,-3 ,-2 ,0 ,-3 ,11 ,-2 ,2 ,-3 ,0 ,-1 ,-2 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,0 ,-1 ,-1 ,-1 ,-1 ,-4 ,-2 ,-1 ,-1 ,0 ,0 ,0 ,-1 ,-2 ,-1 ,-1 ,-1 ,-2 ,-3 ,-2 ,-3 ,-2 ,3 ,-3 ,2 ,-1 ,0 ,-2 ,-1 ,-1 ,-2 ,-4 ,-3 ,-1 ,-2 ,-2 ,-2 ,0 ,-1 ,2 ,-1 ,7 ,-2 ,-1 ,1 ,-3 ,1 ,4 ,-3 ,-2 ,0 ,-3 ,0 ,1 ,-3 ,-1 ,0 ,-4 ,-1 ,3 ,0 ,0 ,-1 ,0 ,-2 ,-3 ,-1 ,-2 ,4 };
	char* seq1 = NULL;
	char* seq2 = NULL;

	GetBlosum("blosum62.csv", blosum);
	//PrintBlosum(blosum);

	switch (argc)
	{
	case 3:
		seq1 = ReadFasta(argv[1]);
		seq2 = ReadFasta(argv[2]);

		//GetBlosum("blosum62.csv", blosum);
		//PrintBlosum(blosom);

		//NeedlemanWunsch("IPGAWD", "VGAWAD", -8, blosum);
		//NeedlemanWunsch("AAAAVGAWADDDD", "VGAWAD", -8, blosum);
		NeedlemanWunsch(seq1, seq2, -4, blosum);

		return 0;

	case 5:
		if (!strcmp(argv[1], "-g"))
		{
			seq1 = ReadFasta(argv[3]);
			seq2 = ReadFasta(argv[4]);

			GetBlosum(argv[2], blosum);
			//PrintBlosum(blosom);

			//NeedlemanWunsch("IPGAWD", "VGAWAD", -8, blosum);
			//NeedlemanWunsch("AAAAVGAWADDDD", "VGAWAD", -8, blosum);
			NeedlemanWunsch(seq1, seq2, -4, blosum);

			return 0;
		}
	default:
		printf("\nNeedleman-Wunsch Alignment Algorithm");
		printf("\n************************************");
		printf("\nImplemented by Haoyang Cheng, Sep 30 2017.\n\n\n");
		printf("NWAlignment [-g gradMat] seq1 seq2\n\n");
		printf("  seq1\t\tThe 1st sequence to be aligned. Expect a .fasta containing only one sequence.\n");
		printf("  seq2\t\tThe 2nd sequence to be aligned. Expect a .fasta containing only one sequence.\n");
		printf("  -g gradMat\tSpecify the grading matrix. See blosum62.csv for format.\n            \tAmino acid order determined by NCBI. Using BLOSUM-62 by default.\n\n\n\n");
		break;

	}
}

#pragma region Needleman-Wunsch

void NeedlemanWunsch(char* seq1, char* seq2, int penalty, int* scoreMat)
{
	int matWidth = strlen(seq1) + 1;
	int matHeight = strlen(seq2) + 1;
	int* gradeMatrix = (int*)malloc(sizeof(int) * matWidth * matHeight);
	int* dirMatrix = (int*)malloc(sizeof(int) * matWidth * matHeight);
	int startHeight = 0;
	int startWidth = 0;
	char* resPath = (char*)malloc(sizeof(char)*(strlen(seq1) + strlen(seq2) + 1));
	char* resSeq1 = (char*)malloc(sizeof(char)*(strlen(seq1) + strlen(seq2) + 1));
	char* resSeq2 = (char*)malloc(sizeof(char)*(strlen(seq1) + strlen(seq2) + 1));

	// Initializing the matrices
	for (int i = 0; i < matWidth * matHeight; i++)
	{
		gradeMatrix[i] = 0;
		dirMatrix[i] = 0;
	}

	for (int i = 0; i < strlen(seq1) + strlen(seq2) + 1; i++)
	{
		resPath[i] = '\0';
		resSeq1[i] = '\0';
		resSeq2[i] = '\0';
	}

	// Initializing the first row of matrices
	// In dirMatrix, -1 means from the left, 1 means from the top, 0 means from the diaganol element.

	gradeMatrix[0] = 0;
	dirMatrix[0] = 0;

	for (int i = 1; i < matWidth; i++)
		gradeMatrix[i] = gradeMatrix[i - 1] + penalty;

	for (int i = 1; i < matHeight; i++)
		gradeMatrix[i * matWidth] = gradeMatrix[(i - 1) * matWidth] + penalty;


	// The real work is here.

	CalMatrices(matHeight, matWidth, penalty, gradeMatrix, dirMatrix, scoreMat, seq1, seq2);
	FindStartPoint(gradeMatrix, matWidth, matHeight, &startWidth, &startHeight);

	//PrintMatrix(gradeMatrix, matWidth, matHeight);
	//printf("\n\n\n");
	//PrintMatrix(dirMatrix, matWidth, matHeight);

	FindPath(dirMatrix, startWidth, startHeight, matWidth, matHeight, resPath);
	GetAlignment(resPath, seq1, seq2, resSeq1, resSeq2);

	ReverseString(resSeq1);
	ReverseString(resSeq2);

	// Printing the results.
	//PrintResultResp(seq1, seq2, resSeq1, resSeq2);
	PrintResultComp(seq1, seq2, resSeq1, resSeq2);

	// Being a good guy: just releasing the memory. 

	if (!gradeMatrix)
	{
		free(gradeMatrix);
		gradeMatrix = NULL;
	}
	if (!dirMatrix)
	{
		free(dirMatrix);
		dirMatrix = NULL;
	}
	if (!resPath)
	{
		free(resPath);
		resPath = NULL;
	}
	if (!resSeq1)
	{
		free(resSeq1);
		resSeq1 = NULL;
	}
	if (!resSeq2)
	{
		free(resSeq2);
		resSeq2 = NULL;
	}
}

// Calculating the values of elem's. in grade/dirMatrix
void CalMatrices(int matHeight, int matWidth, int penalty, int* gradeMatrix, int* dirMatrix, int* scoreMat, char* seq1, char* seq2)
{
	for (int j = 1; j < matHeight; j++)
	{
		for (int i = 1; i < matWidth; i++)
		{
			int leftElem = gradeMatrix[j * matWidth + i - 1];
			int topElem = gradeMatrix[(j - 1) * matWidth + i];
			int diagElem = gradeMatrix[(j - 1) * matWidth + i - 1];
			int scoreItem = scoreMat[(seq2[j - 1] - 'A') * BLOSUM_SIZE + (seq1[i - 1] - 'A')];

			if (topElem + penalty > leftElem + penalty)
			{
				if (topElem + penalty > diagElem + scoreItem)
				{
					gradeMatrix[j * matWidth + i] = topElem + penalty;
					dirMatrix[j * matWidth + i] = 1;
				}
				else
				{
					gradeMatrix[j * matWidth + i] = diagElem + scoreItem;
					dirMatrix[j * matWidth + i] = 0;
				}
			}
			else
			{
				if (leftElem + penalty > diagElem + scoreItem)
				{
					gradeMatrix[j * matWidth + i] = leftElem + penalty;
					dirMatrix[j * matWidth + i] = -1;
				}
				else
				{
					gradeMatrix[j * matWidth + i] = diagElem + scoreItem;
					dirMatrix[j * matWidth + i] = 0;
				}
			}
		}
	}
}

// Find where to start to find a alignment in the matrices.
void FindStartPoint(int* matrix, int width, int height, int* startWidth, int* startHeight)
{

	int maxX = width - 1;
	int maxY = height - 1;
	int maxValue = matrix[height * width - 1];

	for (int i = 1; i < height + 1; i++)
	{
		if (matrix[i * width - 1] > maxValue)
		{
			maxX = width - 1;
			maxY = i - 1;
			maxValue = matrix[i * width - 1];
		}
	}

	for (int i = 0; i < width; i++)
	{
		if (matrix[(height - 1) * width + i] > matrix[maxY * width + maxX])
		{
			maxX = i;
			maxY = height - 1;
			maxValue = matrix[(height - 1) * width + i];
		}
	}

	*startHeight = maxY;
	*startWidth = maxX;
}

// After getting the dirMatrix and where to start, find the path across the matrix
void FindPath(int * dirMatrix, int startWidth, int startHeight, int matWidth, int matHeight, char* path)
{
	int i = startWidth;
	int j = startHeight;
	int k = 0;
	int currDir = 0;

	for (int l = 0; l < matWidth - 1 - startWidth; l++)
	{
		path[k] = 'l';
		k++;
	}

	for (int m = 0; m < matHeight - 1 - startHeight; m++)
	{
		path[k] = 'u';
		k++;
	}

	if (startWidth == matWidth - 1 && startHeight == matHeight - 1)
	{
		//path[k] == 'd';
	}

	while (i >= 0 && j >= 0)
	{
		currDir = dirMatrix[j * matWidth + i];
		switch (currDir)
		{
		case 1:
			j--;
			path[k] = 'u';
			break;
		case 0:
			i--;
			j--;
			path[k] = 'd';
			break;
		case -1:
			i--;
			path[k] = 'l';
			break;
		}
		k++;
	}
}

// From path, get the aligned sequences.
void GetAlignment(char* path, const char* seq1, const char* seq2, char* resChar1, char* resChar2)
{
	int seq1Ptr = strlen(seq1) - 1;
	int seq2Ptr = strlen(seq2) - 1;
	int i;

	for (i = 0; i < strlen(path) && seq1Ptr >= 0 && seq2Ptr >= 0; i++)
	{
		switch (path[i])
		{
		case 'u':
			resChar1[i] = '-';
			resChar2[i] = seq2[seq2Ptr--];
			break;
		case 'l':
			resChar2[i] = '-';
			resChar1[i] = seq1[seq1Ptr--];
			break;
		case 'd':
			resChar1[i] = seq1[seq1Ptr--];
			resChar2[i] = seq2[seq2Ptr--];
			break;
		default:
			break;
		};
	}

	while (seq1Ptr >= 0)
	{
		resChar1[i] = seq1[seq1Ptr--];
		resChar2[i++] = '-';
	}

	while (seq2Ptr >= 0)
	{
		resChar2[i] = seq2[seq2Ptr--];
		resChar1[i++] = '-';
	}


	//for (int i = strlen(path) - 1; i >= 0, seq1Ptr < strlen(seq1), seq2Ptr < strlen(seq2); i--)
	//{
	//	switch (path[i])
	//	{
	//	case 'u':

	//	default:
	//		break;
	//	}
	//}
}

void ReverseString(char * theString)
{
	char temp = '\0';
	int length = strlen(theString);

	for (int i = 0; i < length / 2 - 1; i++)
	{
		temp = theString[i];
		theString[i] = theString[length - 1 - i];
		theString[length - 1 - i] = temp;
	}
}

void PrintResultResp(const char * seq1, const char * seq2, const char * resSeq1, const char * resSeq2)
{
	printf("\nNeedleman-Wunsch Alignment");
	printf("\n**************************");
	printf("\n\n ---  Input Sequences  --- \n\n");
	printf("1: %s\n", seq1);
	printf("2: %s\n", seq2);
	printf("\n\n --- Aligned Sequences --- \n\n");
	printf("1: %s\n", resSeq1);
	printf("2: %s\n", resSeq2);
	printf("\n\n");
}

void PrintResultComp(const char * seq1, const char * seq2, const char * resSeq1, const char * resSeq2)
{
	int i = 0;
	int j = 0;

	printf("\n\n\nNeedleman-Wunsch Alignment Algorithm");
	printf("\n************************************");
	printf("\nImplemented by Haoyang Cheng, Sep 30 2017.\n\n");
	printf("\n\n ---  Input Sequences  --- \n\n");

	printf("1: \n%c", seq1[0]);
	for (int i = 1; i < strlen(seq1); i++)
	{
		if (!(i % NEW_LINE_EVERY_CHARS))
			printf("\n");
		printf("%c", seq1[i]);
	}

	printf("\n\n2: \n%c", seq2[0]);
	for (int i = 1; i < strlen(seq2); i++)
	{
		if (!(i % NEW_LINE_EVERY_CHARS))
			printf("\n");
		printf("%c", seq2[i]);
	}

	printf("\n\n\n\n\n --- Aligned Sequences --- \n\n");

	while (i < strlen(resSeq1) || j < strlen(resSeq2))
	{
		if (i < strlen(resSeq1))
		{
			if (i == 0)
			{
				printf("1: %d\t", i + 1);
			}
			printf("%c", resSeq1[i++]);
			if (i % NEW_LINE_EVERY_CHARS == 0 || i >= strlen(resSeq1))
			{
				printf("\n2: %d\t%c", j, resSeq2[j++]);
				while (j > 0 && j < strlen(resSeq2) && j % NEW_LINE_EVERY_CHARS != 0)
				{
					printf("%c", resSeq2[j++]);
				}
				if (i < strlen(resSeq1))
					printf("\n\n1: %d\t%c", i, resSeq1[i++]);
			}
		}
	}
	printf("\n\n\n\n");
}

void PrintMatrix(int* matrix, int width, int height)
{
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			printf("%d\t", matrix[i * width + j]);
		}
		printf("\n");
	}
}

#pragma endregion


#pragma region FileOperations

char* ReadFasta(char* filename)
{
	FILE* fastaFile = NULL;
	int fileLength = -1;
	char* theString;
	char* currLine = (char*)malloc(sizeof(char) * FASTA_LINE_MAX_LENGTH);
	int writePtr = 0;
	int readPtr = 0;

	fopen_s(&fastaFile, filename, "r");
	if (fastaFile)
	{
		fseek(fastaFile, 0, SEEK_END);
		fileLength = ftell(fastaFile) + 1;
		fseek(fastaFile, 0, SEEK_SET);
	}
	else
		return NULL;

	theString = (char*)malloc(sizeof(char) * fileLength);
	for (int i = 0; i < fileLength; i++)
	{
		theString[i] = '\0';
	}

	fread(theString, sizeof(char), fileLength - 1, fastaFile);

	if (theString[0] == '>')
		while (theString[readPtr++] != '\n');
	for (; theString[readPtr] != '\0'; readPtr++)
	{
		if (theString[readPtr] == '\n')
			continue;
		theString[writePtr++] = theString[readPtr];
	}
	for (; writePtr <= readPtr; writePtr++)
	{
		theString[writePtr] = '\0';
	}

	return theString;
}

#pragma endregion


#pragma region GettingBlosumMatrix

void GetBlosum(char * filename, int * resptr)
{
	FILE* fp = NULL;
	char temp[MAX_LINE_LENGTH] = "";
	int initials[MAX_ARRAY_SIZE] = { 0 };
	int elementNum = 0;
	int lineNum = 0;

	fopen_s(&fp, filename, "r");
	if (!fp)
	{
		return;
	}

	while (!feof(fp))
	{
		ClearString(temp, MAX_LINE_LENGTH);
		fgets(temp, MAX_LINE_LENGTH, fp);
		elementNum = FindInitials(temp, ',', MAX_LINE_LENGTH, initials, MAX_ARRAY_SIZE);
		ChangeDelim(temp, ',', MAX_LINE_LENGTH);
		for (int i = 0; i < BLOSUM_SIZE; i++)
		{
			resptr[(aminoArray[lineNum] - 'A') * BLOSUM_SIZE + aminoArray[i] - 'A'] = atoi(temp + initials[i]);
		}
		lineNum++;
	}

	fclose(fp);
}

void ClearString(char* stringPtr, int length)
{
	for (int i = 0; i < length && stringPtr[i] != '\0'; i++)
	{
		stringPtr[i] = '\0';
	}
}

void ChangeDelim(char* stringPtr, char delim, int length) {
	for (int i = 0; i < length && stringPtr[i] != '\0'; i++)
	{
		if (stringPtr[i] == delim)
		{
			stringPtr[i] = '\0';
		}
	}
}

int FindInitials(char* stringPtr, char delim, int length, int* resArray, int maxSize)
{
	int j = 1;
	resArray[0] = 0;

	for (int i = 0; i < length && j < maxSize && stringPtr[i] != '\0'; i++)
	{
		if (stringPtr[i] == delim)
		{
			resArray[j] = i + 1;
			j++;
		}
	}

	return j;
}

void PrintBlosum(int* blosum)
{
	for (int i = 0; i < BLOSUM_SIZE; i++)
	{
		printf("%c\n", i + 'A');
		for (int j = 0; j < BLOSUM_SIZE; j++)
		{
			printf("%d\t", blosum[i * BLOSUM_SIZE + j]);
		}
		printf("\n");
	}
}

#pragma endregion
