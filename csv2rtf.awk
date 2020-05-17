BEGIN {FS=","
    print "{\\rtf1\\ansi"
#    print "\\paperw11909\\paperh16834\\margl1800\\margr1800\\margt1440\\margb1440\\ftnbj\\aenddoc\\ftnrstcont\\aftnrstcont\\ftnnar\\aftnnrlc"
    print "\\sectd"
    printf "{\\trowd\\trqc\n"
}

{
    # We hope the first record contains the right number of columns 
    if (NR==1) { for (i=1;i<=NF;i++) {printf("\\cellx%d\n",2000*i)}; 
    }
    printf("\\pard\\intbl\\plain\\fnil\\fs22 "); 
    for (i=1;i<=NF-1;i++) {
	printf ("%s\\cell\n",$i);	
    }
    printf ("%s\\cell\\row\n",$NF);
}
END{
    printf "}\n}\n"

} 
