/**
 * @file parallel_blat.cpp
 * @brief A parallel version of blat.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.4
 * @date 2011-09-19
 */

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "misc/blat_record.h"
#include "misc/options_description.h"
#include "sequence/sequence_io.h"
#include "misc/utils.h"
#include "sequence/sequence.h"


using namespace std;

OptionsDescription desc;
int num_threads = omp_get_max_threads();
string ref_filename;
string query_filename;
vector<string> split_files;
set<string> valid_contigs;
double similar = 0.95;

const int MaxLine = (1 << 20);
char line[MaxLine];
char buf[MaxLine];

void CreateFile(const string &filename)
{
    string cmd = FormatString("cat /dev/null > %s", filename.c_str());
    system(cmd.c_str());
}

void RemoveFile(const string &filename)
{
    string cmd = FormatString("rm %s", filename.c_str());
    system(cmd.c_str());
}

void Append(const string &in_file, const string &out_file)
{
    string cmd = FormatString("cat %s >> %s", in_file.c_str(), out_file.c_str());
    system(cmd.c_str());
}

void SplitSequenceFile()
{
    FastaWriter writers[num_threads];
    for (int i = 0; i < num_threads; ++i)
        writers[i].Open(split_files[i]);

    FastaReader reader(query_filename);
    int index = 0;
    Sequence seq;
    string comment;
    while (reader.Read(seq, comment))
    {
        if (valid_contigs.find(comment) == valid_contigs.end())
        {
            if (index/num_threads % 2 == 0)
                writers[index % num_threads].Write(seq, comment);
            else
                writers[num_threads-1 - index % num_threads].Write(seq, comment);
            ++index;
        }
    }
}

void RunBlat(const string &option)
{
    vector<pid_t> id(num_threads);
    for (int i = 0; i < num_threads; ++i)
    {
        pid_t pid = fork();
        if (pid == 0)
        {
            string split_file = split_files[i];
            string cmd = FormatString("blat %s %s %s %s.blat > /dev/null", 
                    option.c_str(), ref_filename.c_str(), split_file.c_str(), split_file.c_str());
            system(cmd.c_str());
            exit(0);
        }
        else
            id[i] = pid;
    }

    for (int i = 0; i < num_threads; ++i)
        waitpid(id[i], NULL, 0);
}

void ParallelBlat(const string &option)
{
    SplitSequenceFile();
    RunBlat(option);
    for (int i = 0; i < num_threads; ++i)
    {
        Append(split_files[i] + ".blat", query_filename + ".blat");
        RemoveFile(split_files[i]);
        RemoveFile(split_files[i] + ".blat");
    }

    FILE *fblat = OpenFile(query_filename + ".blat", "rb");
    while (fgets(line, MaxLine, fblat) != NULL)
    {
        BlatRecord record;
        record.Parse(line);

        if (record.match_count > similar * record.query_length
                && record.match_count > similar * abs(record.ref_to - record.ref_from)
            )
        {
            valid_contigs.insert(record.query_name);
        }
    }
    fclose(fblat);
}

int main(int argc, char *argv[])
{
    desc.AddOption("num_threads", "", num_threads, "number of threads");
    desc.AddOption("similar", "", similar, "similarity");

    try
    {
        desc.Parse(argc, argv);
        
        if (argc < 3)
            throw logic_error("not enough parameters");
    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "parallel_blat - use blat to alignment parallely." << endl;
        cerr << "Usage: parallel_blat ref.fa query.fa" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    ref_filename = argv[1];
    query_filename = argv[2];

    split_files.resize(num_threads);
    for (int i = 0; i < num_threads; ++i)
        split_files[i] = FormatString("%s.split%d", query_filename.c_str(), i);
    CreateFile(query_filename + ".blat");

    deque<string> options;
//    options.push_back(" -noHead -tileSize=18 -minMatch=40 -maxGap=0 -maxIntron=1000 -minIdentity=95 -minScore=100 ");
//    options.push_back(" -noHead -tileSize=18 -minMatch=15 -maxGap=0 -maxIntron=1000 -minIdentity=95 -minScore=100 ");
//    options.push_back(" -noHead -tileSize=18 -minMatch=4 ");
    options.push_back(" -noHead ");

    for (unsigned i = 0; i < options.size(); ++i)
        ParallelBlat(options[i]);

    return 0;
}

