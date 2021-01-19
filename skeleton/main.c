#include<stdio.h>
#include<string.h>

#include "db_index.h"
#include "var_aligner.h"


static int usage(void)
{
	fprintf(stderr, "\n");
    fprintf(stderr, "Program:   %s\n", PACKAGE_NAME); 
    fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
    fprintf(stderr, "Contact:   %s\n", CONTACT); 

    fprintf(stderr, "Usage:     %s <Command> [Options]\n", PACKAGE_NAME); 
	fprintf(stderr, "Command: ");
	fprintf(stderr, "\tindex	index reference sequence\n");
	fprintf(stderr, "\t\taln	align read to reference\n");
	fprintf(stderr, "\n");

	return 1;
}

int main(int argc, char *argv[])
{
	int r = 1;
	if (argc < 2)	return usage();
	if (strcmp(argv[1], "index") == 0)	r = build_db_index(argc, argv);
	else if (strcmp(argv[1], "aln") == 0)	r = run_aligner(argc, argv);
	else if (strcmp(argv[1], "--help") == 0)	return usage();
	else {
		fprintf(stderr, "[Waring!!!] wrong command: '%s'\n", argv[1]);
		return 1;
	}
	//if (!r)
		//fprintf(stderr, "[Main:] Real time: %.3f sec; CPU: %.3f sec\n", realtime() - realtime0, cputime());
		//fprintf(stderr, "succeed\n");

	return 0;
}
