#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<pthread.h>
#include<sys/types.h>
#include<stdbool.h>
#define THREAD_MAX 4
#define TestCases 11
int MAX;
int *a,*b; //arrays to be sorted
void print()
{
	for(int i=0;i<MAX;++i)
		printf("%d ",a[i]);
	printf("\n");
}
void Swap(int *a,int *b)
{
	int temp = *a;
	*a = *b;
	*b = temp;
}
void merge(int lo,int mid,int size)
{
	int s1 = mid-lo+1;
	int s2 = size-mid;
	int a1[s1],a2[s2];
	for(int i=0;i<s1;++i) a1[i] = a[lo+i];
	for(int i=0;i<s2;++i) a2[i] = a[mid+i+1];
	int i=0,j=0,k=lo;
	while(i<s1 && j<s2)
	{
		if(a1[i]>a2[j]) a[k++] = a2[j++];
		else a[k++] = a1[i++];
	}
	while(i < s1) a[k++] = a1[i++];
	while(j < s2) a[k++] = a2[j++];
}
void merge_sort(int lo,int hi)
{
	int m = (lo+hi)/2;
	if(lo < hi)
	{
	
	merge_sort(lo,m);
	merge_sort(m+1,hi);
	merge(lo,m,hi);
	}
}
void *t_merge_sort(void * args)
{
	int part = (int*)args;
	int lo = part * (MAX/4);
	int hi = (part+1) * (MAX/4)-1;
	int m = (lo+hi)/2;
	if(lo < hi)
	{
	merge_sort(lo,m);
	merge_sort(m+1,hi);
	merge(lo,m,hi);
	}
}
int partition(int low, int high)
{
    int pivot = a[high];
    int i = (low - 1);
 
    for (int j = low; j <= high - 1; j++) {
 
        if (a[j] < pivot) {
            i++;
            Swap(&a[i], &a[j]);
        }
    }
    Swap(&a[i + 1], &a[high]);
    return (i + 1);
}
void quickSort(int low, int high)
{
    if (low < high) {
        int pi = partition(low, high);
        quickSort(low, pi - 1);
        quickSort(pi + 1, high);
    }
}
void *t_quick_sort(void *args)
{
	int part = (int*)args;
	int low = part * (MAX/THREAD_MAX);
	int high = (part+1)*(MAX/THREAD_MAX)-1;
	if(low < high)
	{
		int pi = partition(low, high);
        quickSort(low, pi - 1);
        quickSort(pi + 1, high);
	}
}
int findMax(int low,int high)
{
	int rvalue = a[low];
	for(int i=low;i<high;++i) rvalue = rvalue>a[i]?rvalue:a[i];
	return rvalue;
}
void CountSortAscending(int low,int high,int exp)
{
	int size = high-low;
	int output[size];
	int count[10] = {0};
	for(int i=low;i<high;++i)
    {
		count[a[i]/exp%10]++;
    }
	for(int i=1;i<10;++i)
    {
		count[i] += count[i-1];
    }
	for(int i=high-1;i>=low;--i)
    {
		output[--count[a[i]/exp%10]] = a[i];
    }
	for(int i=0;i<size;++i)
    {
		a[low+i] = output[i];
    }
}
void RadixSort(int low,int high)
{
	int max = findMax(low,high);
	for(int exponent = 1;max/exponent > 0;exponent *= 10)
	{
		CountSortAscending(low,high,exponent);
	}
}
void *t_RadixSort(void *args)
{
	int part = (int*)args;
	int low = part * (MAX/THREAD_MAX);
	int high = (part+1)*(MAX/THREAD_MAX);
	RadixSort(low,high);
}
void insertion_sort()
{
    int i, key, j;
    for (i = 1; i < MAX; i++) {
        key = a[i];
        j = i - 1;
        while (j >= 0 && a[j] > key) {
            a[j + 1] = a[j];
            j--;
        }
        a[j + 1] = key;
    }
}
void *t_insertion_sort(void *args)
{
    int i=(int*)args;
    int start=i*(MAX/THREAD_MAX),
    	end=start+(MAX/THREAD_MAX)-1,key,j;
    	i=start;
    while(i<=end)
    {
    	key=a[i];
    	j=i-1;
        while (j >= start && a[j] > key) {
            a[j + 1] = a[j];
            j--;
        }
        a[j + 1] = key;
        i++;
    }
}
int getGap(int i)
{
	i = (i*10)/13;
	if(i<1)
	return 1;
	return i;
}
void CombSort(int low,int high)
{
	int gap = (high-low);
	bool swap = true;
	while(gap != 1 || swap == true)
	{
		swap = false;
		gap = getGap(gap);
		for(int i=low;i<high-gap;++i)
		{
			if(a[i]>a[gap+i])
			{
				Swap(&a[i],&a[gap+i]);
				swap = true;
			}
		}
	}
}
void ShellSort(int low,int high)
{
	int size = high-low;
	for(int gap = high/2;gap>0;gap /= 2)
	{
		for(int i=gap;i<high;++i)
		{
			int key = a[i];
			int j;
			for(j = i;j>=gap && a[j-gap] > key;j-=gap)
			{
				a[j] = a[j-gap];
			}
			a[j] = key;
		}
	}
}
void *t_ShellSort(void *args)
{
	int part = (int*)args;
	int low = part * (MAX/THREAD_MAX);
	int high = (part+1)*(MAX/THREAD_MAX);
	ShellSort(low,high);
}
void *t_CombSort(void *args)
{
	int part = (int*)args;
	int low = part * (MAX/THREAD_MAX);
	int high = (part+1)*(MAX/THREAD_MAX);
	CombSort(low,high);
}
void bubble_sort()
{
    int i, j;
    for (i = 0; i < MAX - 1; i++)
        for (j = 0; j < MAX - i - 1; j++)
        	if (a[j] > a[j + 1])
        		Swap(&a[j], &a[j + 1]);
}
void *t_bubble_sort(void *args)
{
    int i=(int*)args;
    int start=i*(MAX/THREAD_MAX),
    	end=start+(MAX/THREAD_MAX)-1;
    for (int i = start; i < end; i++)
       for(int j=0;j<end-i;j++)
       		if (a[j+start] > a[j+start+ 1])
        		Swap(&a[j+start], &a[j+start+1]);
}
void merge_threads(){
	merge(0,(MAX/2-1)/2,MAX/2-1);
	merge(MAX/2,(MAX/2+MAX-1)/2,MAX-1);
	merge(0,(MAX-1)/2,MAX-1);
}
void make_array()
{
	for(int i=0;i<MAX;++i) b[i] = rand()%100;
}
void set_array()
{
	for(int i=0;i<MAX;++i) a[i] = b[i];
}
int main()
{
	pthread_t p[THREAD_MAX];
	//arrays to store time taken for normal sorting and multithreaded sorting
	float nm[TestCases], tm[TestCases]; //merge sort
	float nq[TestCases], tq[TestCases]; //quick sort
	float ni[TestCases], ti[TestCases]; //insertion sort
	float nr[TestCases], tr[TestCases]; //Radix Sort
	float nc[TestCases], tc[TestCases]; //Comb Sort
	float ns[TestCases], ts[TestCases]; //Shell Sort
	float nb[TestCases], tb[TestCases]; //bubble sort
	float overhead[7][TestCases]; 
	clock_t t1, t2, t3;
	srand(time(NULL));
	MAX=32;
	for(int count=1;count<=TestCases;++count)
	{
		a = (int*)malloc(MAX*sizeof(int));
		b = (int*)malloc(MAX*sizeof(int));
		make_array();
		
		//Merge Sort using process
		set_array();
		t1 = clock();
		merge_sort(0,MAX-1);
		t2 = clock();
		nm[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		
		//Quick Sort using process
		set_array();
		t1 = clock();
		quickSort(0,MAX-1);
		t2 = clock();
		nq[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		
		//Insertion Sort using process
		set_array();
		t1 = clock();
		insertion_sort();
		t2 = clock();
		ni[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		
		//RadixSort using process
		set_array();
		t1 = clock();
		RadixSort(0,MAX);
		t2 = clock();
		nr[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;

		//CombSort using process
		set_array();
		t1 = clock();
		CombSort(0,MAX);
		t2 = clock();
		nc[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		
		//ShellSort using process
		set_array();
		t1 = clock();
		ShellSort(0,MAX);
		t2 = clock();
		ns[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		
		//Bubble Sort using process
		set_array();
		t1 = clock();
		bubble_sort();
		t2 = clock();
		nb[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		
		//Merge Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<4;++i) pthread_create(&p[i],0,t_merge_sort,(void*)i); 
		for(int i=0;i<4;++i) pthread_join(p[i],NULL);
		t2 = clock();
		merge_threads(); //merge the 4 subarrays into one
		t3 = clock();
		tm[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		overhead[0][count-1] = (t3-t2)/(double)CLOCKS_PER_SEC;
		
		//Quick Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<THREAD_MAX;++i) pthread_create(&p[i],0,t_quick_sort,(void*)i);
		for(int i=0;i<THREAD_MAX;++i) pthread_join(p[i],NULL);
		t2 = clock();
		merge_threads(); //merge the 4 subarrays into one
		t3 = clock();
		tq[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		overhead[1][count-1] = (t3-t2)/(double)CLOCKS_PER_SEC;
	
		//Insertion Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<4;++i) pthread_create(&p[i],0,t_insertion_sort,(void*)i);
		for(int i=0;i<4;++i) pthread_join(p[i],NULL);
		t2 = clock();
		merge_threads(); //merge the 4 subarrays into one
		t3 = clock();
		ti[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		overhead[2][count-1] = (t3-t2)/(double)CLOCKS_PER_SEC;
		
		//Radix Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<4;++i) pthread_create(&p[i],0,t_RadixSort,(void*)i);
		for(int i=0;i<4;++i) pthread_join(p[i],NULL);
		t2 = clock();
		merge_threads(); //merge the 4 subarrays into one
		t3 = clock();
		tr[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		overhead[3][count-1] = (t3-t2)/(double)CLOCKS_PER_SEC;
		
		//Comb Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<4;++i) pthread_create(&p[i],0,t_CombSort,(void*)i);
		for(int i=0;i<4;++i) pthread_join(p[i],NULL);
		t2 = clock();
		merge_threads(); //merge the 4 subarrays into one
		t3 = clock();
		tc[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		overhead[4][count-1] = (t3-t2)/(double)CLOCKS_PER_SEC;
		
		//Shell Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<4;++i) pthread_create(&p[i],0,t_ShellSort,(void*)i);
		for(int i=0;i<4;++i) pthread_join(p[i],NULL);
		t2 = clock();
		merge_threads(); //merge the 4 subarrays into one
		t3 = clock();
		ts[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		overhead[5][count-1] = (t3-t2)/(double)CLOCKS_PER_SEC;
		
		//Bubble Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<4;++i) pthread_create(&p[i],0,t_bubble_sort,(void*)i);
		for(int i=0;i<4;++i) pthread_join(p[i],NULL);
		t2 = clock();
		merge_threads(); //merge the 4 subarrays into one
		t3 = clock();
		tb[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		overhead[6][count-1] = (t3-t2)/(double)CLOCKS_PER_SEC;
		
		free(a);
		free(b);
		MAX = MAX*2;
	}
	printf("\n");
	FILE *ptr;
	ptr = fopen("MergeSortResults.txt","w");
	int index=5;
	fprintf(ptr,"Size,Process,Threads\n");
	printf("----------------MERGE SORT----------------------\n");
	printf("SIZE\tPROCESS\t\tTHREADS\t\tMERGING OVERHEAD\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\t\t%f\n", index+i,nm[i],tm[i],overhead[0][i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,nm[i],tm[i]);
	}
	fclose(ptr);
	ptr = fopen("QuickSortResults.txt","w");
	index=5;
	fprintf(ptr,"Size,Process,Threads\n");
	printf("\n----------------QUICK SORT----------------------\n");
	printf("SIZE\tPROCESS\t\tTHREADS\t\tMERGING OVERHEAD\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\t\t%f\n", index+i,nq[i],tq[i],overhead[1][i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,nq[i],tq[i]);
	}
	fclose(ptr);
	ptr = fopen("InsertionSortResults.txt","w");
	index=5;
	fprintf(ptr,"Size,Process,Threads\n");
	printf("\n----------------INSERTION SORT----------------------\n");
	printf("SIZE\tPROCESS\t\tTHREADS\t\tMERGING OVERHEAD\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\t\t%f\n", index+i,ni[i],ti[i],overhead[2][i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,ni[i],ti[i]);
	}
	fclose(ptr);
	ptr = fopen("RadixSortResults.txt","w");
	index=5;
	fprintf(ptr,"Size,Process,Threads\n");
	printf("\n----------------RADIX SORT----------------------\n");
	printf("SIZE\tPROCESS\t\tTHREADS\t\tMERGING OVERHEAD\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\t\t%f\n", index+i,nr[i],tr[i],overhead[3][i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,nr[i],tr[i]);
	}
	fclose(ptr);
	ptr = fopen("CombSortResults.txt","w");
	index=5;
	fprintf(ptr,"Size,Process,Threads\n");
	printf("\n----------------COMB SORT----------------------\n");
	printf("SIZE\tPROCESS\t\tTHREADS\t\tMERGING OVERHEAD\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\t\t%f\n", index+i,nc[i],tc[i],overhead[4][i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,nc[i],tc[i]);
	}
	fclose(ptr);
	ptr = fopen("ShellSortResults.txt","w");
	index=5;
	fprintf(ptr,"Size,Process,Threads\n");
	printf("\n----------------SHELL SORT----------------------\n");
	printf("SIZE\tPROCESS\t\tTHREADS\t\tMERGING OVERHEAD\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\t\t%f\n", index+i,ns[i],ts[i],overhead[5][i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,ns[i],ts[i]);
	}
	fclose(ptr);
	ptr = fopen("BubbleSortResults.txt","w");
	index=5;
	fprintf(ptr,"Size,Process,Threads\n");
	printf("\n----------------BUBBLE SORT----------------------\n");
	printf("SIZE\tPROCESS\t\tTHREADS\t\tMERGING OVERHEAD\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\t\t%f\n", index+i,nb[i],tb[i],overhead[6][i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,nb[i],tb[i]);
	}
	fclose(ptr);
	return 0;
}

