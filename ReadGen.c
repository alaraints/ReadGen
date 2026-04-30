#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <regex.h> 
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>

    
    int readLength = 100;   // -l Read Length  
    int numReads = 1000000; // -n number of reads to generate
    double p = 0;           // -m mutation probability % 
    char** fnamesOut=NULL;  // -o outfiles
    int paired = 0;         // -p paired output
    char* OutPath = NULL;   // -O where the outfiles go (Optional. By default the outfiles go to the same folder as input files.)
    char* title = NULL;     // -t title tag for title lines
    int readlimit = 0;      // -L 0 reads the whole files
    char** fnamesInput=NULL;// -i input file(s). Many files can be entered, gap separated. Fasta format only
    int ForwardOnly = 0;    // -F only forward reads output
    int ReverseOnly = 0;    // -R only reverse reads output
    int BothWays = 1;       // -B output reads both ways mixed  
    int InputCount=0;
    int outputCount=0;
    char* fname1 = NULL;  
    char* outfileR=NULL; 
    char* outfile2=NULL;

    char* msg = "\nReadGen fastq read generator, \n© Alar Aints 2026\n\n Usage:\n\n -L [int] readLimit\n -O outPath\n -i <File(s).fasta>\n -n n reads\n -l reads length>\n -t title_tag>\n -m mutation probability\n -p paired output\n -F/R/B reads out orientation\n";

void parse_args(int argc, char** argv) {
    int opt;
    while ((opt = getopt(argc, argv, "hl:t:L:O:n:m:po:i:FRB")) != -1) {
        switch (opt) {
        case 'h':
            printf("%s\n",msg);
            exit(EXIT_SUCCESS); 
        case 'i': {
            fnamesInput = calloc(argc, sizeof(char*));
            fnamesInput[0]= optarg;
            InputCount++;
            // Collect following arguments that are not options
            while (optind < argc && argv[optind][0] != '-') {
                fnamesInput[InputCount] = argv[optind];
                InputCount++;
                optind++;
            }
            break;
        }
        case 'L':
            readlimit = atoi(optarg);
            break; 
        case 'F':
            ForwardOnly = 1;
            ReverseOnly = 0;
            BothWays = 0;
            break;
        case 'R':
            ReverseOnly = 1;
            ForwardOnly = 0;
            BothWays = 0;
            break;
        case 'B':
            BothWays = 1;
            ForwardOnly = 0;
            ReverseOnly = 0;
            break;
        case 'n':
            numReads = atoi(optarg);
            break;
        case 'l':
            readLength = atoi(optarg);
            break; 
        case 't':
            title = optarg;
            break;  
        case 'm':
            p = atof(optarg);
            break; 
        case 'p':
            paired = 1;
            break;   
        case 'O':
            OutPath = optarg;
            break; 
        case 'o': {
            fnamesOut = calloc(argc, sizeof(char*));
            fnamesOut[0]= optarg;
            outputCount++;
            // Collect following arguments that are not options
            while (optind < argc && argv[optind][0] != '-') {
                fnamesOut[outputCount] = argv[optind];
                outputCount++;
                optind++;
            }
            break;  
        }           
        default:
            fprintf(stderr, "Usage: ReadGen\n%s\n",msg);
            exit(EXIT_FAILURE);
        }
    }
}

const char PT[] = "ACGTN";
char* revTrans(char* seq) {
    int seq_len = strlen(seq);
    char* rt = malloc((seq_len + 1) * sizeof(char));
    int i = seq_len - 1;
    for (; i >= 0; i--) {
        char r = seq[i];
        int j = 0;
        for (; j < 5; j++) {
            if (PT[j] == r) {
                break;
            }
        }
        rt[seq_len - 1 - i] = (j == 4) ? r : PT[3 - j];
    }
    rt[seq_len] = '\0';
    return rt;
}

char mutate_base(char b) {
    char newb = PT[rand() % 5];
    return newb;
}

uint64_t filestatRes[3];

uint64_t* fileStat(char* file) {
    
    uint64_t count = 0;
    char c;
    int cc = 0; // character count per line
    //int ll = 0; // line length
    //int maxll = 0;
    uint64_t fl = 0; // fasta line count
    uint64_t tc = 0; // total characters
    size_t line_length = 0;
    ssize_t read;
    char* line = NULL;
    // Open the file
    FILE *fp = fopen(file, "r");

    // Check if file is open
    if (fp == NULL) {
        fprintf(stderr,"Could not open file %s\n", file);
        return 0;
        exit(1);
    }

    // Count the number of lines
    while ((read = getline(&line, &line_length, fp)) != -1 ) {
        if (line[0] == '>') {
            fl++;
            count++;
            } // fasta count
        else{
            cc = strlen(line)-1;
            tc+=cc;
            count++;
        }
    }
    filestatRes[0]=tc;          // total characters
    filestatRes[1]=count;       // total lines
    filestatRes[2]=fl;          // fasta lines
    // Close the file
    fclose(fp);
    // fprintf(stderr,"Number of lines: %d", count);
    // fprintf(stderr,"Longest line: %d bp\n",maxll);
    return filestatRes;
}

char* openFileInBuffer ( char* fileIn){
    FILE* file = fopen(fileIn, "r");
    if (file == NULL) {
        fprintf(stderr,"Error opening file\n");
        
    }
         // --- 1. Find file size ---
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    rewind(file);  // Go back to start

        // --- 2. Allocate buffer ---
    char *buffer = malloc(size + 2);
    if (!buffer) {
        perror("Memory allocation failed");
        fclose(file);
        
    }
         // --- 3. Read entire file into memory ---
    size_t read_bytes = fread(buffer, 1, size, file);
    fclose(file);
    buffer[read_bytes] = '\0';  // Null-terminate to be safe
        // --- 4. Process line-by-line in memory ---
    return buffer;
}

char* openFasta ( char* fileIn){

    FILE* file = fopen(fileIn, "r");
    if (file == NULL) {
        fprintf(stderr,"Error opening file\n");
        return 0;
    }
         // --- 1. Find file size ---
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    rewind(file);  // Go back to start

        // --- 2. Allocate buffer ---
    char *buffer = malloc(size + 3);
    if (!buffer) {
        perror("Memory allocation failed");
        fclose(file);
        return 0;
    }
         // --- 3. Read entire file into memory ---
    size_t line_length;
    ssize_t read;
    char* line = NULL;
    size_t pos = 0;
    while ((read = getline(&line, &line_length, file)) != -1) {
        if (line[0] == '>'){buffer[pos++] ='\n';
            for (int i = 0; line[i] != '\0'; i++) {
                buffer[pos++] = line[i];
            }
            //fprintf(stderr,"Line %s",line);
        }
        else { 
            for (int i = 0; line[i] != '\0'; i++) {
                if (isalpha(line[i])){
                    if(line[i] >= 'a' && line[i] <= 'z'){line[i] -= 32;} // ensure upper case
                    buffer[pos++] = line[i];
                }
            }
            //fprintf(stderr,"Line %s",line);
        }        
    }
    fclose(file);
    free(line);
    buffer[pos++] = '\n'; 
    buffer[pos++] = '\0';  // Null-terminate to be safe
    return buffer;
}

char* deriveTitleLine (char* fname1){
    size_t lLen = strlen(fname1);
    title = calloc(lLen+1,sizeof(char));
    int dotLen = -1;
    for (int a = 0; a < lLen; a++) {
        int i = lLen - a - 1;
        if (fname1[i] == '.') {
            dotLen = a;
            break;
        }   
    }
    int lastSlashIndex = -1;
        for (int i = lLen; i >= 0; i--) {
            if (fname1[i] == '/') {
                lastSlashIndex = i;
                break;
            }   
        }       
    strncpy(title,fname1+lastSlashIndex+1,lLen-lastSlashIndex-1-dotLen-1);
    size_t tLen = strlen(title);
    for(int i=0;i<tLen;i++){if(title[i]==' '){title[i]='_';}}
    return title;
}

int testFileType (char* fname1){
    FILE *filep = fopen(fname1, "r");
    int firstchar;
    int fileType = 1;
    do {
        firstchar = fgetc(filep);
    } while (firstchar == '\n' || firstchar == '\r' || firstchar == ' ' || firstchar == '\t');
    if (firstchar == EOF) {
        printf("Empty file or error\n");
    } else if (firstchar == '>') {
        fileType=0;
        //disableReset = 1;
        //cycle_step = 1;
        //printf("FASTA format (starts with '>')\n");
    } else if (firstchar == '@') {
        fileType=1;
        // printf("FASTQ format (starts with '@')\n");
    } else {
        printf("Unknown format\n");
    }
    fclose(filep);
    return fileType;
}
int dir_exists_and_writable(const char *path) {
    struct stat st;

    // Check if it exists
    if (stat(path, &st) != 0) {
        return 0;  // does not exist
    }

    // Check if it's a directory
    if (!S_ISDIR(st.st_mode)) {
        return 0;  // not a directory
    }

    // Check if writable
    if (access(path, W_OK) != 0) {
        return 0;  // not writable
    }

    return 1;  // exists, is directory, and writable
}

size_t fasta_root_len(const char* fname) {
    char* fasta = strstr(fname, ".fasta");
    if (!fasta) { fasta = strstr(fname, ".fna"); }
    if (!fasta) { fasta = strstr(fname, ".fa"); }
    if (!fasta) { fasta = strstr(fname, ".txt"); }
    if (fasta) { return (size_t)(fasta - fname); }
    return strlen(fname);
}

int main (int argc, char **argv) {

    parse_args (argc, argv);
    fprintf(stderr, "Running: %s\n",argv[0]);
    if(outputCount == 2){
        paired=1;
        outfileR = fnamesOut[0];
        outfile2 = fnamesOut[1];
    }
    size_t outLen=0; 
    if(OutPath != NULL && OutPath[outLen-1] != '/'){
        OutPath=strdup(OutPath);
        OutPath=realloc(OutPath,(outLen+4)*sizeof(char));
        strncat(OutPath,"/",1);
        outLen=strlen(OutPath);        
    }
    if (OutPath && (! dir_exists_and_writable(OutPath))) {
            mkdir(OutPath, 0777);
    }
    //fprintf(stderr,"outCount: %d\n",outputCount);
    //fprintf(stderr,"paire: %d\n",paired);
    if(outputCount == 1 && paired == 1){  // derive the second filename, check for "R1"
        outfileR = strdup(fnamesOut[0]);
        outfile2 = strdup(fnamesOut[0]);
        //fprintf(stderr,"out1: %s\n",outfileR);
        //    fprintf(stderr,"out2: %s\n",outfile2);
        size_t o2Len = strlen(outfile2);
        char* R1 = strstr(outfile2,"R1");
        if(R1){size_t Rlen = R1 - outfile2;
            outfile2[Rlen+1]='2';
        }
        else{
            char* fq = strstr(fnamesOut[0],"fastq");
            outfileR=calloc((o2Len+9),sizeof(char));
            outfile2=calloc((o2Len+9),sizeof(char));
            if(fq){ 
                size_t rootlen=fq-fnamesOut[0];
                strncpy(outfileR,fnamesOut[0],rootlen);
                strncpy(outfile2,fnamesOut[0],rootlen);
                strncat(outfileR,"R1",2);
                strncat(outfile2,"R2",2);
                strncat(outfileR,".fastq",6);
                strncat(outfile2,".fastq",6);
            }
            else{
                strncpy(outfileR,fnamesOut[0],o2Len);
                strncpy(outfile2,fnamesOut[0],o2Len);
                strncat(outfileR,"R1",2);
                strncat(outfile2,"R2",2);
                strncat(outfileR,".fastq",6);
                strncat(outfile2,".fastq",6);
            }
            //fprintf(stderr,"out1: %s\n",outfileR);
            //fprintf(stderr,"out2: %s\n",outfile2);
        }
        int lastSlashIndex = -1;
        for (int i = o2Len; i >= 0; i--) {
            if (fnamesOut[0][i] == '/') {
                lastSlashIndex = i;
                break;
            }   
        }
        //fprintf(stderr,"out2: %s\n",outfile2);
        if(OutPath){
            char* outfileNew = calloc(outLen + o2Len+9, sizeof(char));
            strncpy(outfileNew,OutPath,outLen);
            strncat(outfileNew,outfileR, o2Len+9);
            outfileR = outfileNew;
            char* outfileNew2 = calloc(outLen + o2Len+9, sizeof(char));
            strncpy(outfileNew2,OutPath,outLen);
            strncat(outfileNew2,outfile2, o2Len+9);
            outfile2 = outfileNew2;

        }
        //fprintf(stderr,"out2: %s\n",outfile2);
    }
    
    if(outputCount == 0 ){
        char* fasta = NULL;
        size_t oRLen = strlen(fnamesInput[0]);
        size_t rootlen = fasta_root_len(fnamesInput[0]);
        if(! OutPath){
                       // derive 1 outfile name                       
            outfileR = calloc((oRLen+22),sizeof(char));
            strncpy(outfileR,fnamesInput[0],rootlen);
            strncat(outfileR,".R1.fastq",10);
            fprintf(stderr,"out1: %s\n",outfileR);
            if (paired==1) {            // derive 2 outfile names
                outfile2 = calloc((oRLen+22),sizeof(char));
                strncpy(outfile2,fnamesInput[0],rootlen);
                strncat(outfile2,".R2.fastq",10);
            }
            
            fprintf(stderr,"out2: %s\n",outfile2);
        } else {
            
            int lastSlashIndex = -1;
            for (int i = oRLen; i >= 0; i--) {
                if (fnamesInput[0][i] == '/') {
                    lastSlashIndex = i;
                    break;
                }   
            }
            size_t tailLen = rootlen - lastSlashIndex +1;
                 // derive 1 outfile name  
                 
            outfileR = calloc(outLen+rootlen + 22,sizeof(char));
            strncpy(outfileR,OutPath,outLen);
            strncat(outfileR,fnamesInput[0] + lastSlashIndex + 1,tailLen );
            strncat(outfileR,".R1.fastq",14);

            if (paired==1) {            // derive 2 outfile names
                outfile2 = calloc(outLen+tailLen + 22,sizeof(char));
                strncpy(outfile2,OutPath,outLen);
                strncat(outfile2,fnamesInput[0] + lastSlashIndex + 1,tailLen );
                strncat(outfile2,".R2.fastq",14);
            }
            
        }

    }
    if(outputCount > 2){ 
        fprintf(stderr,"Max two files can be generated\n");
        exit(1);
    }

    int deriveTitle = 0;
    if(title == NULL){deriveTitle=1;}


  // cycle through input files
  for ( int fi = 0; fi<InputCount; fi+=1){ // start of batch cycle
    fname1=fnamesInput[fi];
    //fprintf(stderr,"Processing: %s \n",fnamesInput[fi]);
    if(deriveTitle){title=deriveTitleLine(fname1);}

    int fileType = testFileType(fname1);
    if(fileType != 0){continue;}         // files not in fasta format will be skipped

  
    if(readlimit == 0){ // get the inputfile linecount via wc -l    
        fileStat(fname1);
        readlimit = (uint32_t)filestatRes[2];
    }
    //fprintf(stderr,"Readlimit: %d\n",readlimit);
    //if(readlimit > 100000){
    //    fprintf(stderr,"Readlimit: %.3f M\n",((float)readlimit)/1000000);
    //}

    char** TITLE = calloc(readlimit+1, sizeof(char*));
    char** seqs = calloc(readlimit+1, sizeof(char*));
    size_t* lengths = calloc(readlimit+1, sizeof(size_t));

    char* FastaBuffer = openFasta(fname1);
    //fprintf(stderr,"Fasta Open\n");
    int m = -1; 
    char* line = FastaBuffer;
    for (char *p = FastaBuffer; *p; ++p) {
        if (*p == '\n') {
            *p = '\0';  // Replace newline with terminator
            size_t lLen = strlen(line);
            //fprintf(stderr,"lLen: %zu\n",lLen);
            if(lLen < 2){line = p + 1;continue; }
            if (line[0] == '>') {
                m+=1; 
                if(m == readlimit){break;}
                TITLE[m]=line;
                //fprintf(stderr,"%d\t%s\n",m,TITLE[m]);
            }
            else{
                seqs[m]=line;
                lengths[m]=lLen;
                //fprintf(stderr,"%s\t%d\t%zu\n",line,m,lengths[m]);
            }
            line = p + 1;
        }
    }
    fprintf(stderr,"Imported %s, %d fasta entries; ",fname1,m+1);
    size_t bp = 0;
    for(int e =0;e<m+1;e++){bp+=lengths[e];}
    fprintf(stderr," %zu bp\n",bp); 
    // end of fasta import 
   
    //fprintf(stderr,"OutCount: %d \n",outputCount);
    char* read = malloc( readLength+1 * sizeof(char));
    char* read2 = NULL;
    
    FILE* fileoutR = fopen(outfileR, "a");
    if (fileoutR == NULL) {
        fprintf(stderr,"Error writing outfile: %s \n",outfileR);
        return 1;
    }
    FILE* fileout2 = NULL;
    if(paired){
        fileout2 = fopen(outfile2, "a");
        if (fileout2 == NULL) {
            fprintf(stderr,"Error writing outfile: %s \n",outfile2);
            return 1;
        }
        read2 = malloc( readLength+1 * sizeof(char));
    }

    for(uint64_t r = 0; r < numReads; r++) {
        
        // pick random sequence
        int s = 0;
        if( m > 0){
            s = rand() % m;
        }
        if(lengths[s] <= 3* readLength ){r--; continue;}
        int start = rand() % (lengths[s] - readLength );
        if(paired ){ start = readLength + (rand() % (lengths[s] - 3*readLength )); }
        memcpy(read, seqs[s] + start, readLength);

        // mutate
        if(p){
            for(int i = 0; i < readLength; i++) {
                double u = (double)rand() / RAND_MAX;
                if(100 * u < p) {
                    read[i] = mutate_base(read[i]);
                }
            }
        }
        read[readLength] = '\0';
        int rev = 0;
        // random strand
        if(ReverseOnly || (BothWays && (rand() % 2)) ){
            read = revTrans(read);
            rev=1;
        }
        int start2 = 0;
        if(paired){ 
            if(rev==0){ start2 = start+ 50 +  (rand() % (300)); if(start > lengths[s]-readLength-1){start2=lengths[s]-readLength-1;}}
            else {start2 = start - 50 - (rand() % (300));if(start2 < 0){start2=0;}}
            memcpy(read2, seqs[s] + start2, readLength); 
            if(p){
                for(int i = 0; i < readLength; i++) {
                    double u = (double)rand() / RAND_MAX;
                    if(100 * u < p) {
                        read2[i] = mutate_base(read2[i]);
                    }
                }
            }
            read2[readLength] = '\0';
            if(!rev){read2=revTrans(read2);}
        }

        // FASTQ output (dummy quality)
        //memcpy(title, TITLE[s] + 1, 42);
        fprintf(fileoutR, "@%s_read%llu/1\n%s\n+\n",title,(unsigned long long)r+1, read);
        for(int i = 0; i < readLength; i++) {fputc('I', fileoutR);}
        fputc('\n', fileoutR);

        if(paired){
            fprintf(fileout2, "@%s_read%llu/2\n%s\n+\n",title,(unsigned long long)r+1, read2);
            for(int i = 0; i < readLength; i++) {fputc('I', fileout2);}
            fputc('\n', fileout2);
        }
        // free(read);read = NULL;

    }
    fclose(fileoutR);
    fclose(fileout2);
    fprintf(stderr,"Outfile: %s \n",outfileR);
    if(paired){fprintf(stderr,"Outfile2: %s \n",outfile2); }
    free(FastaBuffer);FastaBuffer=NULL;
    free(read);read = NULL;
    free(TITLE);TITLE=NULL;
    free(seqs);seqs=NULL;
    free(lengths);lengths=NULL;
    //free(title);title=NULL;
  }
    return 0;

}
