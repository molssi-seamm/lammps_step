/* mdi_nodelay.c -- MolSSI SEAMM: set TCP_NODELAY on MDI's TCP connections.

   MDI writes each message as a header and a body and sets no TCP_NODELAY, so the
   body waits for the header's acknowledgement, which the other side delays by up to
   ~40 ms (Nagle's algorithm and delayed ACKs). For a fast engine that is a large part
   of every step. lammps_step compiles this file in the job's directory and preloads it
   (LD_PRELOAD) into LAMMPS and the engine when they talk over TCP. It wraps connect()
   and accept() and sets TCP_NODELAY on each TCP socket they return; nothing else
   changes. It is not needed once MDI sets TCP_NODELAY itself, and does no harm then.

   Set SEAMM_NODELAY_VERBOSE to see each socket on stderr.
*/
#define _GNU_SOURCE
#include <dlfcn.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>

static void nodelay(int fd)
{
    int one = 1, type = 0;
    socklen_t len = sizeof(type);
    if (fd >= 0 && getsockopt(fd, SOL_SOCKET, SO_TYPE, &type, &len) == 0
        && type == SOCK_STREAM
        && setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one)) == 0
        && getenv("SEAMM_NODELAY_VERBOSE"))
        fprintf(stderr, "mdi_nodelay: TCP_NODELAY on fd %d\n", fd);
}

int connect(int fd, const struct sockaddr *addr, socklen_t len)
{
    static int (*real)(int, const struct sockaddr *, socklen_t);
    if (!real)
        real = dlsym(RTLD_NEXT, "connect");
    int r = real(fd, addr, len);
    if (r == 0 && addr
        && (addr->sa_family == AF_INET || addr->sa_family == AF_INET6))
        nodelay(fd);
    return r;
}

int accept(int fd, struct sockaddr *addr, socklen_t *len)
{
    static int (*real)(int, struct sockaddr *, socklen_t *);
    if (!real)
        real = dlsym(RTLD_NEXT, "accept");
    int r = real(fd, addr, len);
    if (r >= 0)
        nodelay(r); /* MDI passes addr = NULL; nodelay() skips non-TCP sockets */
    return r;
}

int accept4(int fd, struct sockaddr *addr, socklen_t *len, int flags)
{
    static int (*real)(int, struct sockaddr *, socklen_t *, int);
    if (!real)
        real = dlsym(RTLD_NEXT, "accept4");
    int r = real(fd, addr, len, flags);
    if (r >= 0)
        nodelay(r);
    return r;
}
