#ifndef CPEDS_STRING_H
#define CPEDS_STRING_H

#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

inline std::vector<std::string> cpedsSplit(const std::string& text,
                                           const std::string& delimiter) {
    if (delimiter.empty())
        throw std::invalid_argument("String delimiter must not be empty");
    std::vector<std::string> result;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t end = text.find(delimiter, begin);
        result.push_back(
            text.substr(begin, end == std::string::npos ? end : end - begin));
        if (end == std::string::npos)
            break;
        begin = end + delimiter.size();
    }
    return result;
}

inline std::vector<std::string> cpedsSplit(const std::string& text,
                                           char delimiter) {
    return cpedsSplit(text, std::string(1, delimiter));
}

inline std::vector<std::string> cpedsSplitWhitespace(const std::string& text) {
    std::istringstream input(text);
    std::vector<std::string> result;
    std::string token;
    while (input >> token)
        result.push_back(token);
    return result;
}

inline std::string cpedsSimplifyWhitespace(const std::string& text) {
    const std::vector<std::string> words = cpedsSplitWhitespace(text);
    std::string result;
    for (const std::string& word : words) {
        if (!result.empty())
            result += ' ';
        result += word;
    }
    return result;
}

inline std::string cpedsReplaceAll(std::string text, const std::string& search,
                                   const std::string& replacement) {
    if (search.empty())
        throw std::invalid_argument("Search string must not be empty");
    std::size_t position = 0;
    while ((position = text.find(search, position)) != std::string::npos) {
        text.replace(position, search.size(), replacement);
        position += replacement.size();
    }
    return text;
}

inline double cpedsToDouble(const std::string& text, bool* ok = nullptr) {
    char* end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    const bool valid = end != text.c_str() && end != nullptr && *end == '\0';
    if (ok)
        *ok = valid;
    else if (!valid)
        throw std::invalid_argument("Invalid floating-point value: " + text);
    return valid ? value : 0.0;
}

inline long cpedsToLong(const std::string& text) {
    std::size_t consumed = 0;
    const long value = std::stol(text, &consumed);
    if (consumed != text.size())
        throw std::invalid_argument("Invalid integer value: " + text);
    return value;
}

inline int cpedsToInt(const std::string& text) {
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != text.size())
        throw std::invalid_argument("Invalid integer value: " + text);
    return value;
}

#endif
