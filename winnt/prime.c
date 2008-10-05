
#include "windows.h"
#include "main.h"
#include "prime95.h"
#include <direct.h>
#include <math.h>
#include <ctype.h>
#include <dos.h>
#include <fcntl.h>
#include <io.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/timeb.h>

#define PORT	5
#include "gwnum.h"
#include "gwutil.h"
#include "commona.c"
#include "commonb.c"
#include "commonc.c"
#include "ecm.c"
#include "comm95b.c"
#include "comm95c.c"
#include "primenet.c"

void title (char *msg)
{
}

void flashWindowAndBeep ()
{
	MessageBeep (0xFFFFFFFF);
}

/* Return TRUE if we should stop calculating */

int escapeCheck ()
{
	return (THREAD_STOP);
}

void OutputStr (char *buf)
{
	if (DEBUGGING) printf ("%s", buf);
}

void ChangeIcon (
	int	icon_id)
{
}

void BlinkIcon (
	int	duration)
{
}

void BroadcastMessage (
	char	*message)
{
	char	filename[33];
	int	fd;

/* Generate broadcast message file name */

        strcpy (filename, "bcastmsg");
        strcat (filename, EXTENSION);

/* If this is a call to check if a broadcast message exists, then do so */

	if (message == NULL) return;

/* Otherwise, this is a new message - write it to the file */

	fd = _open (filename, _O_TEXT | _O_RDWR | _O_CREAT | _O_APPEND, 0666);
	if (fd < 0) return;
	_write (fd, message, strlen (message));
	_close (fd);
}
