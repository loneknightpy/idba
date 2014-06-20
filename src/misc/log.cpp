/**
 * @file log.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.8
 * @date 2011-12-23
 */

#include "misc/log.h"

#include <dirent.h>
#include <fcntl.h>
#include <pthread.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <stdint.h>
#include <string>

#include "misc/utils.h"

using namespace std;

static const uint32_t kMaxLine = (1 << 20);

LogProcess::LogProcess(const string &log_file)
{
    int pipe_fds[2];
    //int ret = pipe2(pipe_fds, 0);
    int ret = pipe(pipe_fds);
    if (ret == -1)
    {
        throw runtime_error("cannot create pipe");
    }

    pid = fork();
    if (pid == 0)
    {
        close(STDIN_FILENO);
        dup2(pipe_fds[0], STDIN_FILENO);
        close(pipe_fds[0]);

        char line[kMaxLine];
        FILE *flog = OpenFile(log_file, "wb");
        while (fgets(line, kMaxLine, stdin) != NULL)
        {
            fprintf(stdout, "%s", line);
            fprintf(flog, "%s", line);
            fflush(NULL);
        }
        
        exit(0);
    }
    else
    {
        close(STDOUT_FILENO);
        dup2(pipe_fds[1], STDOUT_FILENO);
        close(pipe_fds[1]);
    }
}

LogProcess::~LogProcess()
{
    if (pid != 0)
    {
        int ret = kill(pid, SIGKILL);
        if (ret == -1)
        {
            throw runtime_error("cannot kill log process");
        }

        waitpid(pid, NULL, 0);
    }
}

LogThread::LogThread(const string &log_file)
{
    log_file_ = log_file;

    int pipe_fds[2];
    int ret = pipe(pipe_fds);
    if (ret == -1)
    {
        throw runtime_error("cannot create pipe");
    }

    out_fd = dup(STDOUT_FILENO);
    close(STDOUT_FILENO);
    dup2(pipe_fds[1], STDOUT_FILENO);
    close(pipe_fds[1]);
    in_fd = pipe_fds[0];

    pthread_create(&thread_, NULL, LogThread::LogThreadFunc, (void *)this);
}

LogThread::~LogThread()
{
    pthread_cancel(thread_);
    //close(in_fd);
    pthread_join(thread_, NULL);
}

void *LogThread::LogThreadFunc(void *p)
{
    LogThread &log_thread = *((LogThread *)p);

    FILE *flog = OpenFile(log_thread.log_file_, "wb");

    char line[kMaxLine];
    int len = 0;
    while ((len = read(log_thread.in_fd, line, kMaxLine)) > 0)
    {
        write(log_thread.out_fd, line, len);
        fwrite(line, 1, len, flog);
        fflush(flog);
    }
    fclose(flog);

    return p;
}

