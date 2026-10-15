/*!
 * \file DateTimeData.cpp
 *
 *  Created on: Feb 9, 2018
 *      Author: blew
 */
#include <MscsDateTimeData.h>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace {
std::string datetimeFormat(const std::string& format) {
    std::string result;
    for (std::size_t i = 0; i < format.size();) {
        if (format.compare(i, 4, "yyyy") == 0) {
            result += "%Y";
            i += 4;
        } else if (format.compare(i, 2, "MM") == 0) {
            result += "%m";
            i += 2;
        } else if (format.compare(i, 2, "dd") == 0) {
            result += "%d";
            i += 2;
        } else if (format.compare(i, 2, "hh") == 0) {
            result += "%H";
            i += 2;
        } else if (format.compare(i, 2, "mm") == 0) {
            result += "%M";
            i += 2;
        } else if (format.compare(i, 2, "ss") == 0) {
            result += "%S";
            i += 2;
        } else
            result += format[i++];
    }
    return result;
}

MscsDateTime parseDateTime(const std::string& text, const std::string& format) {
    MscsDateTime result;
    std::istringstream input(text);
    input >> std::get_time(&result.value, datetimeFormat(format).c_str());
    if (input.fail())
        throw std::invalid_argument("Invalid datetime: " + text);
    return result;
}

std::string formatDateTime(const MscsDateTime& datetime,
                           const std::string& format) {
    std::ostringstream output;
    output << std::put_time(&datetime.value, datetimeFormat(format).c_str());
    return output.str();
}

std::vector<std::string> splitWhitespace(const std::string& text) {
    std::istringstream input(text);
    std::vector<std::string> result;
    std::string token;
    while (input >> token)
        result.push_back(token);
    return result;
}
} // namespace

MscsDateTimeData::MscsDateTimeData() {
    //	_f=NULL;
}

MscsDateTimeData::~MscsDateTimeData() {
    //	if (_f!=NULL) { fclose(_f); _f=NULL; }
}

void MscsDateTimeData::load(std::string fname, long dtCol, std::string dtFmt) {
    //	cout << "loading\n";
    ifstream is;
    is.open(fname.c_str(), ios::in);

    if (!is.is_open())
        return;

    _DTdata.dtCol = dtCol;
    _DTdata.dtFmt = dtFmt.c_str();
    _DTdata.fname = fname.c_str();

    long dtSt, dtEn;

    int dtFmt_spaces =
        static_cast<int>(std::count(dtFmt.begin(), dtFmt.end(), ' '));
    dtSt = dtCol;
    dtEn = dtSt + dtFmt_spaces;
    //	cout << "qdtFmt: " << qdtFmt << "\n";
    //	printf("dtSt: %li dtEn: %li\n",dtSt,dtEn);

    double entry;
    std::string line;
    while (getline(is, line)) {
        if (!line.empty() && line[0] != '#') {
            std::vector<std::string> qsl = splitWhitespace(line);

            row_t row;
            std::string dt;
            for (long i = 0; i < qsl.size(); i++) {
                if (i < dtSt)
                    row.data.append(std::stod(qsl[i]));
                else {
                    if (i >= dtSt and i <= dtEn) {
                        dt += qsl[i] + " ";
                    } else {
                        row.data.append(std::stod(qsl[i]));
                    }
                }
            }
            //		cout << dt << "\n";
            //		cout << "dt: " << dt << " qdtFmt: "<< qdtFmt << "\n";
            row.dt = parseDateTime(dt, dtFmt);
            //		cout << "test: " << row.dt.date().year() << " " <<
            //row.dt.date().month() << " " << row.dt.date().day() << " "<<
            //row.dt.time().hour() << ":" << row.dt.time().minute() << ":" <<
            //row.dt.time().second() <<  "\n";
            _DTdata.rows.append(row);
        }
    }
    is.close();

    //	fclose(_f);

    //	cout << "loading done\n";
}

void MscsDateTimeData::printData() {
    for (long i = 0; i < _DTdata.rows.size(); ++i) {
        cout << formatDateTime(_DTdata.rows[i].dt, "yyyy-MM-dd hh:mm:ss")
             << " ";
        for (long j = 0; j < _DTdata.rows[i].data.size(); j++) {
            cout << _DTdata.rows[i].data[j] << " ";
        }
        cout << "\n";
    }
}
/* ********************************************************************************************
 */
double MscsDateTimeData::getFirstValueAfter(std::string dt, long col,
                                            long startFrom, long* idx) {
    MscsDateTime qdt = parseDateTime(dt, _DTdata.dtFmt);
    return getFirstValueAfter(qdt, col, startFrom, idx);
}
/* ********************************************************************************************
 */
double MscsDateTimeData::getFirstValueAfter(MscsDateTime qdt, long col,
                                            long startFrom, long* idx) {
    long i = startFrom;

    while (i < linesCount() and _DTdata.rows[i].dt < qdt) {
        i++;
    }
    double val = 0;
    if (i < linesCount())
        val = _DTdata.rows[i].data[col];
    if (idx != NULL)
        *idx = i;

    return val;
}
/* ********************************************************************************************
 */
double MscsDateTimeData::getFirstValueAfter(double JD, long col, long startFrom,
                                            long* idx) {
    std::string dt =
        cpeds_JDToYMDhms(JD, "%li-%02li-%02li %02li:%02li:%02.0lf");
    MscsDateTime qdt = parseDateTime(dt, "yyyy-MM-dd hh:mm:ss");
    return getFirstValueAfter(qdt, col, startFrom, idx);
}
/* ********************************************************************************************
 */
double MscsDateTimeData::getTimeJD(long idx) {
    const std::tm& time = _DTdata.rows[idx].dt.value;
    double hour = double(time.tm_hour) + double(time.tm_min) / 60 +
                  double(time.tm_sec) / 3600;
    return cpeds_julian_time(time.tm_year + 1900, time.tm_mon + 1, time.tm_mday,
                             hour);
}
