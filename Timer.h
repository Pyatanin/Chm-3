#pragma once

#include <chrono>
#include <iostream>
class Timer
{
public:
	Timer()
	{
		start_timepoint_ = std::chrono::high_resolution_clock::now();
	}
	~Timer()
	{
		stop();
	}

	void stop()
	{
		auto end_timepoint = std::chrono::high_resolution_clock::now();

		auto start =
			std::chrono::time_point_cast<std::chrono::microseconds>(
				start_timepoint_
				)
			.time_since_epoch()
			.count();
		auto end =
			std::chrono::time_point_cast<std::chrono::microseconds>(
				end_timepoint
				)
			.time_since_epoch()
			.count();

		auto duration = end - start;
		double ms = duration * 0.001;

		std::cout << "done for (" << ms << "ms)\n\n";
	}
private:
	std::chrono::time_point<std::chrono::high_resolution_clock> start_timepoint_;
};