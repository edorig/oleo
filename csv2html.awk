#!/usr/bin/awk -f  
BEGIN {
    FS=","
    print "<table>" 
}

{printf("<tr><td>"); 
    for (i=1;i<=NF-1;i++) {
	printf ("%s </td><td> ",$i);	
    }
    printf ("%s </td></tr>\n",$NF);
}
END {print"</table>"}
