#pragma once
#include <iostream>
#include <string>

enum LOG_LEVELS {
	ERROR = 1,
	WARN = 2,
	INFO = 3,
	TRACE = 4,
	DEBUG = 5
};

class Logger {
	private:
		int log_level;
		std::ostream *output_stream;

	public:
		Logger(int level, std::ostream &stream);

		void set_log_level(int level);
		void set_stream(std::ostream &stream);
		void log_msg(std::string msg, int msg_level);
};
