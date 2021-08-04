#pragma once
#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

using namespace std;

class ThreadPool 
{
	private:
		std::vector<std::thread> _threads_;
		std::queue<std::function<void()>> _queue_t_;
		std::mutex _queue_m_;
		std::condition_variable _mes_to_thr_;
		std::atomic<bool> _stop_;

	public:

		template <class Func, class... Args>
			auto exec(Func func, Args... args)->std::future<decltype(func(args...))>
			{
				using type = decltype(func(args...));
				auto task = std::make_shared<std::packaged_task<type()>>(std::packaged_task<type()>(std::move(std::bind(func, args...))));
				std::future<type> res = task->get_future();
				{
					std::unique_lock<std::mutex> lock(_queue_m_);
					if (_stop_)
					{
						throw std::runtime_error("thread terminated");
					}

					_queue_t_.push([task]() { (*task)(); });
				}
				_mes_to_thr_.notify_one();
				return res;
			}

		explicit ThreadPool(size_t size) 
		{
			_stop_ = false;
			for (size_t i = 0; i < size; ++i)
				_threads_.emplace_back([this]() 
						{
						std::function<void()> task;
						while (true) 
						{
						{
						std::unique_lock<std::mutex> lock(_queue_m_);
						_mes_to_thr_.wait(lock, [this] { return !_queue_t_.empty() || _stop_; });
						if (_stop_ && _queue_t_.empty())
						return;
						task = std::move(_queue_t_.front());
						_queue_t_.pop();
						}
						task();
						}
						});
		}

		~ThreadPool() 
		{
			_stop_ = true;
			_mes_to_thr_.notify_all();
			for (std::thread &thread : _threads_)
			{
				thread.join();
			}
		}


		ThreadPool() = delete;
		ThreadPool& operator=(ThreadPool&) = delete;

};
