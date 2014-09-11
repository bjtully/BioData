Getting Started
===============
This will eventually be built into a script, but to start, you will need a file list of the NBCI directory you will be pulling down.
You can do this by copying this script into your favorite text editor and running it with some edits.

```
from ftplib import FTP
ftp = FTP('ftp.ncbi.nih.gov')
ftp.login('anonymous', <EMAIL>)
ftp.cwd("/refseq/release/bacteria")
files = ftp.dir()
print files
```

Be sure to fill in your email address. If you would like complete archaeal genomes, change /refseq/release/bacteria to /refseq/release/archaea.