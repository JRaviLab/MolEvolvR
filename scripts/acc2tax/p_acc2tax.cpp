#include <string>
#include <sstream>
#include <pthread.h>
#include <fstream>
#include <iostream>
#include <time.h>
// #include <unordered_map>
#include <unordered_set>
// #include <set>
// #include <vector>

// 346 seconds: why does it take longer?

using namespace std;

pthread_mutex_t accessions_mutex;
string acc2tax_fp, outfile;

unordered_set<string> accessions;

long total_left = 0;

void *Worker(void *arg);
bool IsPosInt(std::string s);
int main(int argc, char **argv)
{
    // Parse command line arguments
    int NThreads = 1;
    string logS = "./data/acc2tax.log";

    for (int i = 1; i < argc; ++i)
    {
        string cmd = string(argv[i]);
        if (cmd == "-t")
        {
            if ((i + 1) == argc || !IsPosInt(string(argv[i + 1])))
            {
                cout << "Error: -t argument requires a positive integer" << endl;
                return -1;
            }
            else
            {
                string nt = string(argv[++i]);
                unsigned int nti = stoi(nt);
                if (nti == 0)
                {
                    cout << "Error: -t argument requires a positive integer" << endl;
                    return -1;
                }
                else
                {
                    NThreads = nti;
                }
            }
        }
        else if (cmd == "-f")
        {
            if (++i == argc)
            {
                cout << "Error: -f argument requires filepath to accession - taxid mapping" << endl;
                return -1;
            }
            else
            {
                acc2tax_fp = string(argv[i]);
            }
        }
        else if (cmd == "-o")
        {
            if (++i == argc)
            {
                cout << "Error: -o argument requires filepath to the out" << endl;
            }
            else
            {
                outfile = string(argv[i]);
            }
        }
        else if (cmd == "-l")
        {
            if (++i == argc)
            {
                cout << "Error: -l argument requires filepath to the logfile" << endl;
            }
            else
            {
                logS = string(argv[i]);
            }
        }
        else
        {

            ifstream accFile(cmd);
            if (!accFile.is_open())
            {
                cout << "Warning: File " << cmd << " does not exist" << endl;
            }

            else
            {
                string a;
                while (!accFile.eof())
                {

                    accFile >> a;
                    accessions.insert(a);
                    // accessions.push_back(a);
                }
            }
        }
    }

    /* 
    ofstream logfile(logS);
    logfile << "Lines\t"
            << "Seconds\t"
            << "Acc2TaxHits" << endl;
            */

    total_left = accessions.size();

    pthread_t threads[NThreads];
    pthread_mutex_init(&accessions_mutex, NULL);

    time_t total_start, total_end;
    time(&total_start);
    ifstream acc2tax_file(acc2tax_fp);

    // Create Threads
    for (int i = 0; i < NThreads; ++i)
    {
        if (pthread_create(&threads[i], NULL, Worker, (void *)i))
        {
            printf("*** Error creating thread ***\n");
            exit(-1);
        }
    }

    // Join Threads
    for (int i = 0; i < NThreads; ++i)
    {
        if (pthread_join(threads[i], NULL))
        {
            printf("*** Error joining thread ***\n");
            exit(-2);
        }
    }

    time(&total_end);
    cout << "Total Time: " << difftime(total_end, total_start) << " Seconds" << endl;

    return 0;
}

void *Worker(void *arg)
{
    long i = (long)arg;
    string n_acc2tax_fp = acc2tax_fp + to_string(i);
    string outpath = outfile + to_string(i);

    ifstream acc2tax(n_acc2tax_fp);
    if (!acc2tax.is_open())
    {
        cout << "Warning: File " << n_acc2tax_fp << " does not exist" << endl;
        pthread_exit(NULL);
    }

    ofstream out(outpath);
    if (!out.is_open())
    {
        cout << "Warning: Error opening " << outpath << endl;
    }

    string acc, tax;
    while (!acc2tax.eof() && total_left)
    {
        acc2tax >> acc >> tax;

        auto loc = accessions.find(acc);
        if (loc != accessions.end())
        {
            // Item read in is in the set of accessions:
            // Write to output
            out << acc << "\t" << tax << endl;

            pthread_mutex_lock(&accessions_mutex);
            --total_left;
            pthread_mutex_unlock(&accessions_mutex);
            // Deleting actually isn't necessary -- we can instead track how many have been obtained
            // Probably faster b/c reduces busy waiting
            // Delete entry from set
            // accessions.erase(loc);
        }
    }
    acc2tax.close();
    out.close();

    pthread_exit(NULL);
}

bool IsPosInt(std::string s)
{
    // If string contains any characters that are not digits
    // This includes '-' such that negative numbers will also return false
    for (char c : s)
    {
        if (!isdigit(c))
        {
            return false;
        }
    }
    return true;
}