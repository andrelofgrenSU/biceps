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
#include <iostream>
#include <string>

/**
 * @enum LOG_LEVELS
 * @brief Enum representing different log levels for the Logger class.
 * 
 * The log levels are used to categorize the importance of log messages.
 */
enum LOG_LEVELS {
    ERROR = 1, /**< Error log level, for critical errors */
    WARN = 2,  /**< Warning log level, for non-critical issues */
    INFO = 3,  /**< Info log level, for informational messages */
    TRACE = 4, /**< Trace log level, for detailed logs */
    DEBUG = 5  /**< Debug log level, for debugging purposes */
};

/**
 * @class Logger
 * @brief A class responsible for logging messages at different log levels.
 * 
 * The Logger class allows for logging messages of varying levels (e.g., ERROR, WARN, INFO, etc.)
 * and can output these messages to different streams (e.g., console, file).
 */
class Logger {
private:
    int level; /**< Current log level */
    std::ostream *output_stream; /**< Output stream for logging messages */

public:
    /**
     * @brief Default constructor.
     * 
     * Initializes the Logger with the default log level (INFO) and the default output stream (std::cout).
     */
    Logger();

    /**
     * @brief Parameterized constructor to initialize Logger with a specific log level and output stream.
     * 
     * @param level The log level to set for the logger (e.g., ERROR, INFO).
     * @param stream The output stream where log messages will be written (e.g., std::cout, std::ofstream).
     */
    Logger(int level, std::ostream &stream);

    /**
     * @brief Getter for the current log level.
     * 
     * @return The current log level of the logger.
     */
    int log_level();

    /**
     * @brief Sets the log level of the logger.
     * 
     * @param value The log level to set (e.g., ERROR, INFO).
     */
    void set_log_level(int value);

    /**
     * @brief Sets the output stream for the logger.
     * 
     * @param stream The output stream where log messages will be written.
     */
    void set_stream(std::ostream &stream);

    /**
     * @brief Logs a message with a specified log level.
     * 
     * @param msg The message to be logged.
     * @param msg_level The log level for the message (e.g., ERROR, INFO).
     */
    void log_msg(std::string msg, int msg_level);
};
