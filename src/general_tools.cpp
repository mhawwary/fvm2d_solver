

#include"general_tools.h"

// Sorting double array
void QuickSort(double*& szArray,double*& ssArray , int nLower, int nUpper)
{
    // Check for non-base case
    if (nLower < nUpper)
    {
        // Split and sort partitions
        int nSplit = Partition(szArray,ssArray, nLower, nUpper);
        QuickSort(szArray,ssArray, nLower, nSplit - 1);
        QuickSort(szArray,ssArray, nSplit + 1, nUpper);
    }
}

// QuickSort partition implementation
int Partition(double*& szArray,double*& ssArray, int nLower, int nUpper)
{
    // Pivot with first element
    int nLeft = nLower + 1;
    double szPivot = szArray[nLower];
    //double ssPivot=ssArray[nLower];
    int nRight = nUpper;

    // Partition array elements
    double szSwap,ssSwap;
    while (nLeft <= nRight)
    {
        // Find item out of place
        while (nLeft <= nRight && szArray[nLeft]<= szPivot)
            nLeft = nLeft + 1;
        while (nLeft <= nRight && szArray[nRight]>szPivot)
            nRight = nRight - 1;

        // Swap values if necessary
        if (nLeft < nRight)
        {
            szSwap = szArray[nLeft];             ssSwap=ssArray[nLeft];
            szArray[nLeft] = szArray[nRight];    ssArray[nLeft]=ssArray[nRight];
            szArray[nRight] = szSwap;            ssArray[nRight] = ssSwap;
            nLeft = nLeft + 1;
            nRight = nRight - 1;
        }
    }

    // Move pivot element
    szSwap = szArray[nLower];
    szArray[nLower] = szArray[nRight];
    szArray[nRight] = szSwap;

    ssSwap = ssArray[nLower];
    ssArray[nLower] = ssArray[nRight];
    ssArray[nRight] = ssSwap;

    return nRight;
}

// Sorting two double arrays and an ID int one
void QuickSort3(double*& mainArray,double*& A1 , int*& A2 , int nLower, int nUpper)
{
    // Check for non-base case
    if (nLower < nUpper)
    {
        // Split and sort partitions
        int nSplit = Partition3(mainArray,A1,A2,nLower, nUpper);
        QuickSort3(mainArray,A1,A2, nLower, nSplit - 1);
        QuickSort3(mainArray,A1,A2, nSplit + 1, nUpper);
    }
}

// QuickSort partition implementation for three arrays
int Partition3(double*& mainArray,double*& A1 , int*& A2, int nLower, int nUpper)
{
    // Pivot with first element
    int nLeft = nLower + 1;
    double mainPivot =  mainArray[nLower];
    int nRight = nUpper;

    // Partition array elements
    double mainSwap,A1Swap;  int A2Swap;
    while (nLeft <= nRight)
    {
        // Find item out of place
        while (nLeft <= nRight && mainArray[nLeft]<= mainPivot)
            nLeft = nLeft + 1;
        while (nLeft <= nRight && mainArray[nRight]>mainPivot)
            nRight = nRight - 1;

        // Swap values if necessary
        if (nLeft < nRight)
        {
            mainSwap = mainArray[nLeft];     mainArray[nLeft] = mainArray[nRight];     mainArray[nRight] = mainSwap;

            A1Swap= A1[nLeft];  A1[nLeft]=A1[nRight];   A1[nRight] = A1Swap;
            A2Swap= A2[nLeft];  A2[nLeft]=A2[nRight];   A2[nRight] = A2Swap;

            nLeft = nLeft + 1;
            nRight = nRight - 1;
        }
    }

    // Move pivot element
    mainSwap = mainArray[nLower];   mainArray[nLower] = mainArray[nRight];    mainArray[nRight] = mainSwap;

    A1Swap = A1[nLower];    A1[nLower] = A1[nRight];    A1[nRight] = A1Swap;
    A2Swap = A2[nLower];    A2[nLower] = A2[nRight];    A2[nRight] = A2Swap;

    return nRight;
}



