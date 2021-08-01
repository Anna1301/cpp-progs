#include<iostream>
class Allocator{
	private:
		char *start_mem;
		size_t current;
		size_t msize;
	public:
		Allocator()
		{
		start_mem = nullptr;
		current = 0;
		msize = 0;
		}
		~Allocator()
		{
			std::free(start_mem);
		}
		void makeAllocator(size_t maxSize)
		{
			
			if(start_mem != nullptr){
			start_mem = static_cast<char*>(std::realloc(start_mem,maxSize));
			}else{
			start_mem = static_cast<char*>(std::malloc(maxSize* sizeof(char)));
			}

			msize = maxSize;
			current = 0;
		}
		char* alloc(size_t size)
		{
			if(current + size > msize){
			return nullptr;
			} else {
			current += size;
			return start_mem + current - size;
			}
			
		}
		void reset()
		{
			current = 0;
		}

		size_t ret_current(){
		return current;
		}

		size_t ret_msize(){
		return msize;
		}

		char* ret_start_mem(){
		return start_mem;
		}
                

};



