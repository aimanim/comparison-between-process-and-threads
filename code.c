#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<pthread.h>
#include<sys/types.h>
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
	float nm[20], tm[20]; //merge sort
	float nq[20], tq[20]; //quick sort
	float ni[20], ti[20]; //insertion sort
	clock_t t1, t2;
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
		
		//Merge Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<4;++i) pthread_create(&p[i],0,t_merge_sort,(void*)i); 
		for(int i=0;i<4;++i) pthread_join(p[i],NULL);
		merge_threads(); //merge the 4 subarrays into one
		t2 = clock();
		tm[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		
		//Quick Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<THREAD_MAX;++i) pthread_create(&p[i],0,t_quick_sort,(void*)i);
		for(int i=0;i<THREAD_MAX;++i) pthread_join(p[i],NULL);
		merge_threads(); //merge the 4 subarrays into one
		t2 = clock();
		tq[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		
		//Insertion Sort using threads
		set_array();
		t1 = clock();
		//divide array in 4 threads
		for(int i=0;i<4;++i) pthread_create(&p[i],0,t_insertion_sort,(void*)i);
		for(int i=0;i<4;++i) pthread_join(p[i],NULL);
		merge_threads(); //merge the 4 subarrays into one
		t2 = clock();
		ti[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		
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
	printf("SIZE\tPROCESS\t\tTHREADS\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\n", index+i,nm[i],tm[i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,nm[i],tm[i]);
	}
	fclose(ptr);
	ptr = fopen("QuickSortResults.txt","w");
	index=5;
	fprintf(ptr,"Size,Process,Threads\n");
	printf("\n----------------QUICK SORT----------------------\n");
	printf("SIZE\tPROCESS\t\tTHREADS\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\n", index+i,nq[i],tq[i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,nq[i],tq[i]);
	}
	fclose(ptr);
	ptr = fopen("InsertionSortResults.txt","w");
	index=5;
	fprintf(ptr,"Size,Process,Threads\n");
	printf("\n----------------INSERTION SORT----------------------\n");
	printf("SIZE\tPROCESS\t\tTHREADS\n");
	for(int i=0;i<TestCases;++i){
		printf("2^%d\t%f\t%f\n", index+i,ni[i],ti[i]);
		fprintf(ptr,"2^%d,%f,%f\n", index+i,ni[i],ti[i]);
	}
	fclose(ptr);
	return 0;
}
