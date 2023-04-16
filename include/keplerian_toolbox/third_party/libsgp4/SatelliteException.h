/*
 * Copyright 2013 Daniel Warner <contact@danrw.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SATELLITEEXCEPTION_H_
#define SATELLITEEXCEPTION_H_

#include <exception>

/**
 * @brief The exception that the SGP4 class throws upon an error.
 */
class SatelliteException : public std::exception
{
public:
    SatelliteException(const char *message) : message_(message)
    {
    }

    virtual ~SatelliteException(void) noexcept
    {
    }

    virtual const char *what() const noexcept override
    {
        return message_.c_str();
    }

private:
    std::string message_;
};

#endif
