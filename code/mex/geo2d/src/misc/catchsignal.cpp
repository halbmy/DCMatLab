void signalhandler();
void catchsignal(int sig);

#include <stdio.h>
#include <signal.h>
#include <stdlib.h>

#ifndef WIN32
#include <unistd.h>

void signalhandler()
{
  signal(SIGHUP,&catchsignal);
  signal(SIGINT,&catchsignal);
  signal(SIGQUIT,&catchsignal);
  signal(SIGTERM,&catchsignal);
  signal(SIGTSTP,&catchsignal);
  signal(SIGCONT,&catchsignal);
  signal(SIGSEGV,&catchsignal);
}

void catchsignal(int sig)
{
  switch (sig)
  {
  case SIGHUP:
    printf("pausing program execution on signal SIGHUP.\n");
    pause();
    break;
  case SIGTSTP:
    printf("pausing program execution on signal SIGTSTP.\n");
    pause();
    break;
  case SIGCONT:
    printf("resuming program execution on signal SIGCONT.\n");
    break;
  case SIGINT:
    printf("terminating program execution on signal SIGINT.\n");
    exit(0);
    break;
  case SIGTERM:
    printf("terminating program execution on signal SIGTERM.\n");
    exit(0);
    break;
  case SIGQUIT:
    printf("terminating program execution on signal SIGQUIT.\n");
    exit(0);
    break;
  case SIGSEGV:
    printf("terminating program execution on signal SIGSEGV.\n");
    exit(0);
    break;
  default:
    printf("terminating program execution on unhandled signal.\n");
    exit(0);
    break;
  }
}

#else

void signalhandler()
{
  signal(SIGINT,&catchsignal);
  signal(SIGTERM,&catchsignal);
  signal(SIGBREAK,&catchsignal);
}

void catchsignal(int sig)
{
    printf("terminating program execution\n");
    exit(0);
}
#endif

