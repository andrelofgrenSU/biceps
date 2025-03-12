#include <logger.hpp>

Logger::Logger(int level, std::ostream &stream) : level(level), output_stream(&stream) {}
Logger::Logger() : level(WARN), output_stream(&std::cout) {}

int Logger::log_level()
{
    return level;
}

void Logger::set_log_level(int value)
{
    level = value;
}

void Logger::set_stream(std::ostream &stream) {
    output_stream = &stream;
}

void Logger::log_msg(std::string msg, int msg_level) {
    if (msg_level > level)
        return;
    if (output_stream)
        *output_stream << msg << std::endl;
}
