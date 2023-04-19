#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <set>
#include <iterator>
#include <queue>
#include <numeric>
#include <cstdlib>
#include <time.h>
#include <chrono>
#include <iomanip>
#include <string.h>
#define ll long long int
#define ull unsigned long long int
#define all(a) (a).begin(), (a).end()
using namespace std;
using namespace std::chrono;

//ifstream cin("input.txt");
//ofstream cout("output.txt");

void printArray(vector<int> A, int size)
{
    for (auto i = 0; i < size; i++)
        cout << A[i] << " ";
    cout << "\n";
}

/* This function takes last element as pivot, places
the pivot element at its correct position in sorted
array, and places all smaller (smaller than pivot)
to left of pivot and all greater elements to right
of pivot */
int partition(vector<int>& arr, int low, int high)
{
    int pivot = arr[high]; // pivot
    int i = (low - 1); // Index of smaller element and indicates
    // the right position of pivot found so far

    for (int j = low; j <= high - 1; j++) {
        // If current element is smaller than the pivot
        if (arr[j] < pivot) {
            i++; // increment index of smaller element
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return (i + 1);
}

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void quickSort(vector<int>& arr, int low, int high)
{
    if (low < high) {
        /* pi is partitioning index, arr[p] is now
        at right place */
        int pi = partition(arr, low, high);

        // Separately sort elements before
        // partition and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

// Merges two subarrays of array[].
// First subarray is arr[begin..mid]
// Second subarray is arr[mid+1..end]
void merge(vector<int>& array, int const left, int const mid,
    int const right)
{
    auto const subArrayOne = mid - left + 1;
    auto const subArrayTwo = right - mid;

    // Create temp arrays
    auto* leftArray = new int[subArrayOne],
        * rightArray = new int[subArrayTwo];

    // Copy data to temp arrays leftArray[] and rightArray[]
    for (auto i = 0; i < subArrayOne; i++)
        leftArray[i] = array[left + i];
    for (auto j = 0; j < subArrayTwo; j++)
        rightArray[j] = array[mid + 1 + j];

    auto indexOfSubArrayOne
        = 0, // Initial index of first sub-array
        indexOfSubArrayTwo
        = 0; // Initial index of second sub-array
    int indexOfMergedArray
        = left; // Initial index of merged array

    // Merge the temp arrays back into array[left..right]
    while (indexOfSubArrayOne < subArrayOne
        && indexOfSubArrayTwo < subArrayTwo) {
        if (leftArray[indexOfSubArrayOne]
            <= rightArray[indexOfSubArrayTwo]) {
            array[indexOfMergedArray]
                = leftArray[indexOfSubArrayOne];
            indexOfSubArrayOne++;
        }
        else {
            array[indexOfMergedArray]
                = rightArray[indexOfSubArrayTwo];
            indexOfSubArrayTwo++;
        }
        indexOfMergedArray++;
    }
    // Copy the remaining elements of
    // left[], if there are any
    while (indexOfSubArrayOne < subArrayOne) {
        array[indexOfMergedArray]
            = leftArray[indexOfSubArrayOne];
        indexOfSubArrayOne++;
        indexOfMergedArray++;
    }
    // Copy the remaining elements of
    // right[], if there are any
    while (indexOfSubArrayTwo < subArrayTwo) {
        array[indexOfMergedArray]
            = rightArray[indexOfSubArrayTwo];
        indexOfSubArrayTwo++;
        indexOfMergedArray++;
    }
    delete[] leftArray;
    delete[] rightArray;
}

// begin is for left index and end is
// right index of the sub-array
// of arr to be sorted */
void mergeSort(vector<int>& array, int const begin, int const end)
{
    if (begin >= end)
        return; // Returns recursively

    auto mid = begin + (end - begin) / 2;
    mergeSort(array, begin, mid);
    mergeSort(array, mid + 1, end);
    merge(array, begin, mid, end);
}

// To heapify a subtree rooted with node i
// which is an index in arr[].
// n is size of heap
void heapify(vector<int>& arr, int N, int i)
{

    // Initialize largest as root
    int largest = i;

    // left = 2*i + 1
    int l = 2 * i + 1;

    // right = 2*i + 2
    int r = 2 * i + 2;

    // If left child is larger than root
    if (l < N && arr[l] > arr[largest])
        largest = l;

    // If right child is larger than largest
    // so far
    if (r < N && arr[r] > arr[largest])
        largest = r;

    // If largest is not root
    if (largest != i) {
        swap(arr[i], arr[largest]);

        // Recursively heapify the affected
        // sub-tree
        heapify(arr, N, largest);
    }
}

// Main function to do heap sort
void heapSort(vector<int>& arr, int N)
{

    // Build heap (rearrange array)
    for (int i = N / 2 - 1; i >= 0; i--)
        heapify(arr, N, i);

    // One by one extract an element
    // from heap
    for (int i = N - 1; i > 0; i--) {

        // Move current root to end
        swap(arr[0], arr[i]);

        // call max heapify on the reduced heap
        heapify(arr, i, 0);
    }
}

// Function to sort an array using
// insertion sort
void insertionSort(vector<int>& arr, int n)
{
    int i, key, j;
    for (i = 1; i < n; i++)
    {
        key = arr[i];
        j = i - 1;

        // Move elements of arr[0..i-1], 
        // that are greater than key, to one
        // position ahead of their
        // current position
        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

const int RUN = 32;

// This function sorts array from left index to
// to right index which is of size atmost RUN
void insertionSort(vector<int>& arr, int left, int right)
{
    for (int i = left + 1; i <= right; i++)
    {
        int temp = arr[i];
        int j = i - 1;
        while (j >= left && arr[j] > temp)
        {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = temp;
    }
}

// Merge function merges the sorted runs
void mergetim(vector<int>& arr, int l, int m, int r)
{

    // Original array is broken in two parts
    // left and right array
    int len1 = m - l + 1, len2 = r - m;
    vector<int> left(len1), right(len2);
    for (int i = 0; i < len1; i++)
        left[i] = arr[l + i];
    for (int i = 0; i < len2; i++)
        right[i] = arr[m + 1 + i];

    int i = 0;
    int j = 0;
    int k = l;

    // After comparing, we
    // merge those two array
    // in larger sub array
    while (i < len1 && j < len2)
    {
        if (left[i] <= right[j])
        {
            arr[k] = left[i];
            i++;
        }
        else
        {
            arr[k] = right[j];
            j++;
        }
        k++;
    }

    // Copy remaining elements of left, if any
    while (i < len1)
    {
        arr[k] = left[i];
        k++;
        i++;
    }

    // Copy remaining element of right, if any
    while (j < len2)
    {
        arr[k] = right[j];
        k++;
        j++;
    }
}

// Iterative Timsort function to sort the
// array[0...n-1] (similar to merge sort)
void timSort(vector<int>& arr, int n)
{

    // Sort individual subarrays of size RUN
    for (int i = 0; i < n; i += RUN)
        insertionSort(arr, i, min((i + RUN - 1),
            (n - 1)));

    // Start merging from size RUN (or 32).
    // It will merge
    // to form size 64, then 128, 256
    // and so on ....
    for (int size = RUN; size < n;
        size = 2 * size)
    {

        // pick starting point of
        // left sub array. We
        // are going to merge
        // arr[left..left+size-1]
        // and arr[left+size, left+2*size-1]
        // After every merge, we
        // increase left by 2*size
        for (int left = 0; left < n;
            left += 2 * size)
        {

            // find ending point of
            // left sub array
            // mid+1 is starting point
            // of right sub array
            int mid = left + size - 1;
            int right = min((left + 2 * size - 1),
                (n - 1));

            // merge sub array arr[left.....mid] &
            // arr[mid+1....right]
            if (mid < right)
                mergetim(arr, left, mid, right);
        }
    }
}

//Swap function
void swap(int* xp, int* yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void selectionSort(vector<int>& arr, int n)
{
    int i, j, min_idx;
    // One by one move boundary of 
    // unsorted subarray 
    for (i = 0; i < n - 1; i++)
    {
        // Find the minimum element in 
        // unsorted array 
        min_idx = i;
        for (j = i + 1; j < n; j++)
        {
            if (arr[j] < arr[min_idx])
                min_idx = j;
        }
        // Swap the found minimum element 
        // with the first element 
        if (min_idx != i)
            swap(&arr[min_idx], &arr[i]);
    }
}

/* function to sort arr using shellSort */
int shellSort(vector<int>& arr, int n)
{
    // Start with a big gap, then reduce the gap
    for (int gap = n / 2; gap > 0; gap /= 2)
    {
        // Do a gapped insertion sort for this gap size.
        // The first gap elements a[0..gap-1] are already in gapped order
        // keep adding one more element until the entire array is
        // gap sorted 
        for (int i = gap; i < n; i += 1)
        {
            // add a[i] to the elements that have been gap sorted
            // save a[i] in temp and make a hole at position i
            int temp = arr[i];

            // shift earlier gap-sorted elements up until the correct 
            // location for a[i] is found
            int j;
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap)
                arr[j] = arr[j - gap];

            //  put temp (the original a[i]) in its correct location
            arr[j] = temp;
        }
    }
    return 0;
}

// A function to implement bubble sort
void bubbleSort(vector<int>& arr, int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++)

        // Last i elements are already 
        // in place
        for (j = 0; j < n - i - 1; j++)
            if (arr[j] > arr[j + 1])
                swap(arr[j], arr[j + 1]);
}

struct Node
{
    int key;
    struct Node* left, * right;
};

// A utility function to create a new BST Node
struct Node* newNode(int item)
{
    struct Node* temp = new Node;
    temp->key = item;
    temp->left = temp->right = NULL;
    return temp;
}

// Stores inorder traversal of the BST
// in arr[]
void storeSorted(Node* root, vector<int>& arr, int& i)
{
    if (root != NULL)
    {
        storeSorted(root->left, arr, i);
        arr[i++] = root->key;
        storeSorted(root->right, arr, i);
    }
}

/* A utility function to insert a new
   Node with given key in BST */
Node* insert(Node* node, int key)
{
    /* If the tree is empty, return a new Node */
    if (node == NULL) return newNode(key);

    /* Otherwise, recur down the tree */
    if (key < node->key)
        node->left = insert(node->left, key);
    else if (key > node->key)
        node->right = insert(node->right, key);

    /* return the (unchanged) Node pointer */
    return node;
}

// This function sorts arr[0..n-1] using Tree Sort
void treeSort(vector<int>& arr, int n)
{
    struct Node* root = NULL;

    // Construct the BST
    root = insert(root, arr[0]);
    for (int i = 1; i < n; i++)
        root = insert(root, arr[i]);

    // Store inorder traversal of the BST
    // in arr[]
    int i = 0;
    storeSorted(root, arr, i);
}

// A utility function to get maximum value in arr[]
int getMax(vector<int>& arr, int n)
{
    int mx = arr[0];
    for (int i = 1; i < n; i++)
        if (arr[i] > mx)
            mx = arr[i];
    return mx;
}

// A function to do counting sort of arr[] according to
// the digit represented by exp.
void countSort(vector<int>& arr, int n, int exp)
{
    vector<int> output(n);
    int i, count[10] = { 0 };

    // Store count of occurrences in count[]
    for (i = 0; i < n; i++)
        count[(arr[i] / exp) % 10]++;

    // Change count[i] so that count[i] now contains actual
    //  position of this digit in output[]
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];

    // Build the output array
    for (i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }

    // Copy the output array to arr[], so that arr[] now
    // contains sorted numbers according to current digit
    for (i = 0; i < n; i++)
        arr[i] = output[i];
}

// The main function to that sorts arr[] of size n using
// Radix Sort
void radixsort(vector<int>& arr, int n)
{
    // Find the maximum number to know number of digits
    int m = getMax(arr, n);

    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    for (int exp = 1; m / exp > 0; exp *= 10)
        countSort(arr, n, exp);
}

/* Merge the sorted ranges [low, mid1), [mid1,mid2)
and [mid2, high) mid1 is first midpoint
index in overall range to merge mid2 is second
midpoint index in overall range to merge*/
void merge3(vector<int>& gArray, int low, int mid1,
    int mid2, int high, vector<int>& destArray)
{
    int i = low, j = mid1, k = mid2, l = low;

    // Choose smaller of the smallest in the three ranges
    while ((i < mid1) && (j < mid2) && (k < high))
    {
        if (gArray[i] < gArray[j])
        {
            if (gArray[i] < gArray[k])
            {
                destArray[l++] = gArray[i++];
            }
            else
            {
                destArray[l++] = gArray[k++];
            }
        }
        else
        {
            if (gArray[j] < gArray[k])
            {
                destArray[l++] = gArray[j++];
            }
            else
            {
                destArray[l++] = gArray[k++];
            }
        }
    }

    // Case where first and second ranges
    // have remaining values
    while ((i < mid1) && (j < mid2))
    {
        if (gArray[i] < gArray[j])
        {
            destArray[l++] = gArray[i++];
        }
        else
        {
            destArray[l++] = gArray[j++];
        }
    }

    // case where second and third ranges
    // have remaining values
    while ((j < mid2) && (k < high))
    {
        if (gArray[j] < gArray[k])
        {
            destArray[l++] = gArray[j++];
        }
        else
        {
            destArray[l++] = gArray[k++];
        }
    }

    // Case where first and third ranges have
    // remaining values
    while ((i < mid1) && (k < high))
    {
        if (gArray[i] < gArray[k])
        {
            destArray[l++] = gArray[i++];
        }
        else
        {
            destArray[l++] = gArray[k++];
        }
    }

    // Copy remaining values from the first range
    while (i < mid1)
        destArray[l++] = gArray[i++];

    // Copy remaining values from the second range
    while (j < mid2)
        destArray[l++] = gArray[j++];

    // Copy remaining values from the third range
    while (k < high)
        destArray[l++] = gArray[k++];
}


/* Performing the merge sort algorithm on the
given array of values in the rangeof indices
[low, high). low is minimum index, high is
maximum index (exclusive) */
void mergeSort3WayRec(vector<int>& gArray, int low,
    int high, vector<int>& destArray)
{
    // If array size is 1 then do nothing
    if (high - low < 2)
        return;

    // Splitting array into 3 parts
    int mid1 = low + ((high - low) / 3);
    int mid2 = low + 2 * ((high - low) / 3) + 1;

    // Sorting 3 arrays recursively
    mergeSort3WayRec(destArray, low, mid1, gArray);
    mergeSort3WayRec(destArray, mid1, mid2, gArray);
    mergeSort3WayRec(destArray, mid2, high, gArray);

    // Merging the sorted arrays
    merge3(destArray, low, mid1, mid2, high, gArray);
}

void mergeSort3Way(vector<int>& gArray, int n)
{
    // if array size is zero return null
    if (n == 0)
        return;

    // creating duplicate of given array
    vector<int> fArray(n);

    // copying elements of given array into
    // duplicate array
    for (int i = 0; i < n; i++)
        fArray[i] = gArray[i];

    // sort function
    mergeSort3WayRec(fArray, 0, n, gArray);

    // copy back elements of duplicate array
    // to given array
    for (int i = 0; i < n; i++)
        gArray[i] = fArray[i];
}

void solve() {
    srand(time(0));
    int n = 10000000, m = 1;
    int sum_quicksort = 0, sum_mergesort = 0, sum_heapsort = 0, sum_insertion = 0, sum_timsort = 0,
        sum_selection = 0, sum_shellsort = 0, sum_bubble = 0, sum_treesort = 0, sum_radixsort = 0,
        sum_stlsort = 0, sum_3waymerge = 0;
    for (int j = 1; j <= m; j++) {
        vector<int> sample(n);
        for (int i = 0; i < n; i++) {
            int x;
            x = rand();
            x <<= 15;
            x ^= rand();
            x %= 1000001;
            sample[i] = x;
        }
        //printArray(sample, n - 1);

        vector<int> v;

        v = sample;
        //printArray(v, n - 1);
        auto start = high_resolution_clock::now();
        quickSort(v, 0, n - 1);
        auto stop = high_resolution_clock::now();
        //printArray(v, n - 1);
        auto duration = duration_cast<microseconds>(stop - start);
        //cout << duration.count() << "\n";
        sum_quicksort += duration.count();

        v = sample;
        //printArray(v, n - 1);
        start = high_resolution_clock::now();
        mergeSort(v, 0, n - 1);
        stop = high_resolution_clock::now();
        //printArray(v, n - 1);
        duration = duration_cast<microseconds>(stop - start);
        //cout << duration.count() << "\n";
        sum_mergesort += duration.count();

        v = sample;
        //printArray(v, n - 1);
        start = high_resolution_clock::now();
        heapSort(v, n);
        stop = high_resolution_clock::now();
        //printArray(v, n - 1);
        duration = duration_cast<microseconds>(stop - start);
        //cout << duration.count() << "\n";
        sum_heapsort += duration.count();

        //v = sample;
        ////printArray(v, n - 1);
        //start = high_resolution_clock::now();
        //insertionSort(v, n);
        //stop = high_resolution_clock::now();
        ////printArray(v, n - 1);
        //duration = duration_cast<microseconds>(stop - start);
        ////cout << duration.count() << "\n";
        //sum_insertion += duration.count();

        v = sample;
        //printArray(v, n - 1);
        start = high_resolution_clock::now();
        timSort(v, n);
        stop = high_resolution_clock::now();
        //printArray(v, n - 1);
        duration = duration_cast<microseconds>(stop - start);
        //cout << duration.count() << "\n";
        sum_timsort += duration.count();

        //v = sample;
        ////printArray(v, n - 1);
        //start = high_resolution_clock::now();
        //selectionSort(v, n);
        //stop = high_resolution_clock::now();
        ////printArray(v, n - 1);
        //duration = duration_cast<microseconds>(stop - start);
        ////cout << duration.count() << "\n";
        //sum_selection += duration.count();

        v = sample;
        //printArray(v, n - 1);
        start = high_resolution_clock::now();
        shellSort(v, n);
        stop = high_resolution_clock::now();
        //printArray(v, n - 1);
        duration = duration_cast<microseconds>(stop - start);
        //cout << duration.count() << "\n";
        sum_shellsort += duration.count();

        //v = sample;
        ////printArray(v, n - 1);
        //start = high_resolution_clock::now();
        //bubbleSort(v, n);
        //stop = high_resolution_clock::now();
        ////printArray(v, n - 1);
        //duration = duration_cast<microseconds>(stop - start);
        ////cout << duration.count() << "\n";
        //sum_bubble += duration.count();

        v = sample;
        //printArray(v, n - 1);
        start = high_resolution_clock::now();
        treeSort(v, n);
        stop = high_resolution_clock::now();
        //printArray(v, n - 1);
        duration = duration_cast<microseconds>(stop - start);
        //cout << duration.count() << "\n";
        sum_treesort += duration.count();

        v = sample;
        //printArray(v, n - 1);
        start = high_resolution_clock::now();
        radixsort(v, n);
        stop = high_resolution_clock::now();
        //printArray(v, n - 1);
        duration = duration_cast<microseconds>(stop - start);
        //cout << duration.count() << "\n";
        sum_radixsort += duration.count();

        v = sample;
        //printArray(v, n - 1);
        start = high_resolution_clock::now();
        sort(all(v));
        stop = high_resolution_clock::now();
        //printArray(v, n - 1);
        duration = duration_cast<microseconds>(stop - start);
        //cout << duration.count() << "\n";
        sum_stlsort += duration.count();

        v = sample;
        //printArray(v, n - 1);
        start = high_resolution_clock::now();
        mergeSort3Way(v, n);
        stop = high_resolution_clock::now();
        //printArray(v, n - 1);
        duration = duration_cast<microseconds>(stop - start);
        //cout << duration.count() << "\n";
        sum_3waymerge += duration.count();
    }
    cout << "Average quicksort: \n" << (1.0 * sum_quicksort / (m * 1000)) << "\n";
    cout << "Average mergesort: \n" << (1.0 * sum_mergesort / (m * 1000)) << "\n";
    cout << "Average heapsort: \n" << (1.0 * sum_heapsort / (m * 1000)) << "\n";
    cout << "Average insertion: \n" << (1.0 * sum_insertion / (m * 1000)) << "\n";
    cout << "Average timsort: \n" << (1.0 * sum_timsort / (m * 1000)) << "\n";
    cout << "Average selection: \n" << (1.0 * sum_selection / (m * 1000)) << "\n";
    cout << "Average shellsort: \n" << (1.0 * sum_shellsort / (m * 1000)) << "\n";
    cout << "Average bubble: \n" << (1.0 * sum_bubble / (m * 1000)) << "\n";
    cout << "Average treesort: \n" << (1.0 * sum_treesort / (m * 1000)) << "\n";
    cout << "Average radixsort: \n" << (1.0 * sum_radixsort / (m * 1000)) << "\n";
    cout << "Average stlsort: \n" << (1.0 * sum_stlsort / (m * 1000)) << "\n";
    cout << "Average 3waymerge: \n" << (1.0 * sum_3waymerge / (m * 1000)) << "\n";
}

int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0); cout.tie(0);

    int t = 1;
    //cin >> t;
    while (t--)
        solve();
}