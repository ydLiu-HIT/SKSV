#ifndef _DB_INDEX_H
#define _DB_INDEX_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<getopt.h>
#include<unistd.h>
#include<zlib.h>
#include "desc.h"

int db_index_usage(void)
{
    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   %s\n", PACKAGE_NAME); 
    fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
    fprintf(stderr, "Contact:   %s\n", CONTACT); 
    fprintf(stderr, "Usage:     %s index [Options] <Reference> <HashIndexDir>\n", PACKAGE_NAME); 
    fprintf(stderr, "<Reference>            Sequence of reference genome, in FASTA format\n");
    fprintf(stderr, "<HashIndexDir>         The directory to store deBGA index\n");
	fprintf(stderr, "\nfor building deBGA index file. You can get more deBGA information from https://github.com/HongzheGuo/deBGA\n");
	fprintf(stderr, "\n");
	return 1;
}

int get_bin_dir(char *bin, char *dir)
{
	char *end = strrchr(bin, '/');
	if (end == NULL)
		return 1;

	bin[end-bin] = '\0';
	strcpy(dir, bin);
	return 0;
}

int deBGA_index(char *dir, char *ref_fa, char *index_route)
{
	char cmd[1024];
	sprintf(cmd, "%sdeBGA index %s %s", dir, ref_fa, index_route);
	fprintf(stderr, "[db_index] Executing deBGA index ...\n");
	if (system(cmd) != 0)
	{
		fprintf(stderr, "\n[tgs_index] Indexing undoing, deBGA index exit abnormally. \n"); 
		exit(1);
	}
	fprintf(stderr, "[db_index] Done!\n");

    return 1;
}

int build_db_index(int argc, char *argv[])
{
	// clock_t t = clock();
	char dir[1024];
	char *ref_fa = 0;
	char *index_route = 0;

	if(optind + 3 != argc)	return db_index_usage();

	if (!get_bin_dir(argv[0], dir))
	{
		strcat(dir, "/");
	}
	ref_fa = strdup(argv[optind+1]);
	index_route = strdup(argv[optind+2]);

	deBGA_index(dir, ref_fa, index_route);

	return 0;
}

#endif
