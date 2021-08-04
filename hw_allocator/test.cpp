#include<iostream>
#include<cassert>
#include"alctr.cpp"
void DefaultWorkTest()
{
	Allocator A;
	A.makeAllocator(100);
	char* result_0 = A.alloc(0);
	char* result_1= A.alloc(50);
	char* result_2 = A.alloc(20);
	char* result_3 = A.alloc(0);
	assert(result_3-result_0 == 70); 
}

void DefaultCreateTest()
{
	Allocator A;
	A.makeAllocator(20);
	assert(A.ret_start_mem() != nullptr && A.ret_current() == 0 && A.ret_msize() == 20 );

}

void CheckResetTest()
{
	Allocator A;
	A.makeAllocator(200);
	char* result_0 = A.alloc(0);
	char* result_1 = A.alloc(60);
	char* result_2 = A.alloc(70);
	A.reset();
	assert(A.alloc(0) == result_0);
}

void OutOfBufferTest()
{
	Allocator A;
	A.makeAllocator(50);
	char* result = A.alloc(54);
	assert(result == nullptr);
}

void UsingMemAfterReset()
{
	Allocator A;
	A.makeAllocator(20);
	char* result_1 = A.alloc(10);
	A.reset();
	char* result_2 = A.alloc(10);
	//char* result_3 = A.alloc(0);
	assert(result_1 == result_2);
}

void FullMemFill()
{
	Allocator A;
	A.makeAllocator(51);
	char* result_0 =  A.alloc(0);
	char* result_1 = A.alloc(30);
	char* result_2 =  A.alloc(20);
	char* result_3 = A.alloc(0);
	assert(result_3 == result_0 + 50);
}
int main(){
	DefaultWorkTest();
	CheckResetTest();
	OutOfBufferTest();
	UsingMemAfterReset();
	FullMemFill();
	DefaultCreateTest();

	std::cout << "success \n";

	return 0;
}
