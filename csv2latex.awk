BEGIN {FS=","}

{for (i=1;i<=NF-1;i++) {
	printf ("%s & ",$i);	
    }
    printf ("%s \\\\\n",$NF);
}
