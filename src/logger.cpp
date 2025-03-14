/*
 * Copyright (C) 2025 André Löfgren
 *
 * This file is part of Biceps.
 *
 * Biceps is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Biceps is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Biceps. If not, see <https://www.gnu.org/licenses/>.
 */
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
