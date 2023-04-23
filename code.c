#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<pthread.h>
#include<sys/types.h>
#define THREAD_MAX 4
#define TestCases 15
int MAX;
int *a,*b;
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
	float nm[20],tm[20],nq[20],tq[20];
	clock_t t1, t2;
	srand(time(NULL));
	MAX=32;
	for(int count=1;count<=TestCases;++count)
	{
		MAX = MAX*2;
		a = (int*)malloc(MAX*sizeof(int));
		b = (int*)malloc(MAX*sizeof(int));
		make_array();
		printf("\n\n----------------------ArraySize: %d----------------------",MAX);
		set_array();
		t1 = clock();
		merge_sort(0,MAX-1);
		t2 = clock();
		nm[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		printf("\nTime taken for normal merge sort: %f",nm[count-1] );
		
		
		set_array();
		t1 = clock();
		quickSort(0,MAX-1);
		t2 = clock();
		nq[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		printf("\nTime taken for normal quick sort: %f",nq[count-1] );
		
		set_array();
		t1 = clock();
		
		for(int i=0;i<4;++i)
			pthread_create(&p[i],0,t_merge_sort,(void*)i);
		for(int i=0;i<4;++i)
			pthread_join(p[i],NULL);
		
		merge(0,(MAX/2-1)/2,MAX/2-1);
		merge(MAX/2,(MAX/2+MAX-1)/2,MAX-1);
		merge(0,(MAX-1)/2,MAX-1);
		
		t2 = clock();
		tm[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		printf("\nTime taken for multithreaded MergeSort: %f\n", tm[count-1]);
		
		set_array();
		t1 = clock();
		for(int i=0;i<THREAD_MAX;++i) pthread_create(&p[i],0,t_quick_sort,(void*)i);
		for(int i=0;i<THREAD_MAX;++i) pthread_join(p[i],NULL);
		
		merge(0,(MAX/2-1)/2,MAX/2-1);
		merge(MAX/2,(MAX/2+MAX-1)/2,MAX-1);
		merge(0,(MAX-1)/2,MAX-1);
		t2 = clock();
		tq[count-1] = (t2-t1)/(double)CLOCKS_PER_SEC;
		printf("Time taken for multithreaded QuickSort: %f\n", tq[count-1]);
		free(a);
		free(b);
		
	}
	printf("\n");
	FILE *ptr;
	ptr = fopen("MergeSortResults.txt","w");
	int index=6;
	fprintf(ptr,"Size,Process,Threads\n");
	for(int i=0;i<TestCases;++i){
		fprintf(ptr,"2^%d,%f,%f\n", index+i,nm[i],tm[i]);
		MAX*=2;
		}
	fclose(ptr);
	ptr = fopen("QuickSortResults.txt","w");
	index=6;
	fprintf(ptr,"Size,Process,Threads\n");
	for(int i=0;i<TestCases;++i){
		fprintf(ptr,"2^%d,%f,%f\n", index+i,nq[i],tq[i]);
		MAX*=2;
		}
	fclose(ptr);
	return 0;
}
