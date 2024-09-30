/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright © 2024 Domenic Franjic
 *
 * This file is part of ReplicationNowcastingMacroVarsWithSDFM.
 *
 * ReplicationNowcastingMacroVarsWithSDFM is free software: you can redistribute
 * it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.

 * ReplicationNowcastingMacroVarsWithSDFM is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ReplicationNowcastingMacroVarsWithSDFM. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#ifndef EMPIRICS_SETUP
#define EMPIRICS_SETUP

// Including external libraries
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <Eigen>
#ifdef _WIN32
#include <direct.h>
#define mkdir _mkdir
#else
#include <unistd.h>
#endif

namespace EmpiricsSetup {

    /* Create a pseudo date classe for aligning actual data */
    class Date
    {
    private:

        Eigen::MatrixXd dateMatrix;

    public:

        /* Constructor */
        Date(const Eigen::MatrixXd& matrix) { dateMatrix = matrix; }

        /* Function to get the formatted date */
        std::string giveDate(size_t rowIndex) const
        {
            if (rowIndex < dateMatrix.rows())
            {
                std::ostringstream oss;
                oss << std::setfill('0') << std::setw(4) << static_cast<int>(dateMatrix(rowIndex, 0)) << "-"
                    << std::setw(2) << static_cast<int>(dateMatrix(rowIndex, 1)) << "-"
                    << std::setw(2) << static_cast<int>(dateMatrix(rowIndex, 2));
                return oss.str();

            }
            else
            {
                std::cerr << "Index out of bounds." << std::endl;
                return "";
            }
        }

        /* Function to get the column index based on a given date string */
        int giveIndex(const std::string& dateString)
        {
            for (int i = 0; i < dateMatrix.rows(); ++i)
            {
                if (giveDate(i) == dateString)
                {
                    return i;
                }
            }
            return -1;
        }

        /* Convert the date to the respective matrix index */
        Eigen::VectorXi numericDate()
        {
            return Eigen::VectorXi::LinSpaced(dateMatrix.rows(), 0, dateMatrix.rows() - 1);
        }

        /* Create a vector of indices corresponding to dates*/
        Eigen::VectorXi dateSequenceNumeric(const std::string& first_date, const std::string& last_date, const int frequency)
        {
            int start = giveIndex(first_date), end = giveIndex(last_date);
            if (frequency == 4 && start != -1 && end != -1)
            {
                if ((end - start) % 3 != 0)
                {
                    std::cout << "\n\n Error! The difference between start and end does not correspond to a quaterly time series!" << endl;
                    return Eigen::VectorXi::Constant(1, -1);

                }
                Eigen::VectorXi v = Eigen::VectorXi::Zero((end - start) / 3 + 1);
                int ind = 0;
                for (int i = start; i <= end; i += 3)
                {
                    v(ind) = i;
                    ++ind;
                }
                return v;
            }
            else if (frequency == 12 && start != -1 && end != -1)
            {
                std::cout << "diff = " << end - start << '\n';
                return Eigen::VectorXi::LinSpaced(end - start + 1, start, end);
            }
            return Eigen::VectorXi::Constant(1, -1);
        }

        /* Function to get a sequence of dates corresponding to the indices provided by its argument */
        std::vector<std::string> dateSequence(const Eigen::VectorXi& v)
        {
            std::vector<std::string> date_dates;
            for (int i = 0; i < v.size(); ++i)
            {
                int index = v(i);
                if (index >= 0 && index < dateMatrix.rows())
                {
                    date_dates.push_back(giveDate(index));
                }
                else
                {
                    std::cerr << "Error: Index " << index << " is out of bounds or invalid." << std::endl;

                }
            }
            return date_dates;
        }
    };

    /* Create a class to create, save, and laod the parameters of a pseudo-nowcasting exercise */
    class NowcastSetup {
    public:
        int G1; // Grid size for the l1 penalty
        int G2; // Grid size for the l2 penalty
        int T_CV; // Number of observations used int the validation set
        int H; // Numebr of nowcasts computed in the validation set
        int quarterly; // Whether or not the VOI
        std::vector<std::string> fc_dates; // Dates as string

        // Constructor
        NowcastSetup() : G1(0), G2(0), T_CV(0), H(0), quarterly(0) {}

        /* Save function */
        void saveSetup(const std::string& filename) {
            std::ofstream file(filename, std::ios::binary);
            if (file.is_open()) {
                file.write(reinterpret_cast<char*>(&G1), sizeof(G1));
                file.write(reinterpret_cast<char*>(&G2), sizeof(G2));
                file.write(reinterpret_cast<char*>(&T_CV), sizeof(T_CV));
                file.write(reinterpret_cast<char*>(&H), sizeof(H));
                file.write(reinterpret_cast<char*>(&quarterly), sizeof(quarterly));
                int fcSize = fc_dates.size();
                file.write(reinterpret_cast<char*>(&fcSize), sizeof(fcSize));
                for (const auto& str : fc_dates) {
                    int len = str.size();
                    file.write(reinterpret_cast<char*>(&len), sizeof(len));
                    file.write(str.data(), len);
                }
                file.close();
                std::cout << "Data saved to file: " << filename << std::endl;
            }
            else {
                std::cerr << "Unable to open file: " << filename << std::endl;
            }
        }

        /* Load function */
        void loadSetup(const std::string& filename) {
            std::ifstream file(filename, std::ios::binary);
            if (file.is_open()) {
                file.read(reinterpret_cast<char*>(&G1), sizeof(G1));
                file.read(reinterpret_cast<char*>(&G2), sizeof(G2));
                file.read(reinterpret_cast<char*>(&T_CV), sizeof(T_CV));
                file.read(reinterpret_cast<char*>(&H), sizeof(H));
                file.read(reinterpret_cast<char*>(&quarterly), sizeof(quarterly));
                int fcSize;
                file.read(reinterpret_cast<char*>(&fcSize), sizeof(fcSize));
                fc_dates.clear();
                for (int i = 0; i < fcSize; ++i) {
                    int len;
                    file.read(reinterpret_cast<char*>(&len), sizeof(len));
                    std::string str(len, '\0');
                    file.read(&str[0], len);
                    fc_dates.push_back(str);
                }
                file.close();
                std::cout << "Data loaded from file: " << filename << std::endl;
            }
            else {
                std::cerr << "Unable to open file: " << filename << std::endl;
            }
        }

        /* Print function */
        void printInfo(const Date& dates) {
            std::cout << "\n\nThe current standard nowcasting setup for this data set is as follows:\n";
            std::cout << "Size of l1 penalty grid: " << G1 << std::endl;
            std::cout << "Size of l2 penalty grid: " << G2 << std::endl;
            std::cout << "Date of the last CV nowcast: " << dates.giveDate(T_CV) << std::endl;
            std::cout << "Number of CV nowcasts: " << H << std::endl;
            std::cout << "Frequency of the VOI: " << quarterly << std::endl;
            std::cout << "\nForecasting dates: ";
            for (const auto& d : fc_dates) {
                std::cout << d << ", ";
            }
            std::cout << std::endl;
        }
    };

};
#endif /* defined(EMPIRICS_SETUP) */
