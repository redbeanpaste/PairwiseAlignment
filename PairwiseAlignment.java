package pairwisealignment;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;


public class PairwiseAlignment 
{
    public static int count = 1;
    
    void align(FileReader seq, FileReader scoring) throws IOException
    {
        Object[] obj = read(seq, scoring);
         
        int[] scoreModel = (int[]) obj[2];
        StringBuilder seq1 = (StringBuilder) obj[0];
        StringBuilder seq2 = (StringBuilder) obj[1];
        System.out.println("seq1: "+seq1);
        System.out.println("seq2: "+seq2);
        int match = scoreModel[0];
        int mismatch = scoreModel[1];
        int gap = scoreModel[2];
        int n = seq1.length(); //col
        int m = seq2.length(); // row
        int[][] score = new int[m+1][n+1];
        String[][] traceback = new String[m+1][n+1];
        int i = 0;
        int j = 0;
        // inialize the scoring and traceback table
        for(; i <= m; i++) 
        {
            score[i][j] = (i)*gap;
            traceback[i][j] = "U";
        }
        i = 0;
        for(; j <= n; j++) 
        {
            score[i][j] = (j)*gap;
            traceback[i][j] = "L";
        }
        // start calculate score of two sequence
        for (i = 1; i <= m; i++)
        {
            j = 1;
            for(; j <= n; j++)
            {
                if (seq1.charAt(j-1) == seq2.charAt(i-1))
                // if symbol at index i of 1st seq is match with symbol at index j of 2nd seq
                // compare with 1) symbol at index i of 1st seq align with gap(insert to 2nd sequence)
                //              2) symbol at index j of 2nd seq align with gap(insert to 1st sequence)
                // and add traceback symbol at index [i][j]  
                    score[i][j] = max(score[i-1][j]+gap
                                    ,score[i][j-1]+gap
                                    ,score[i-1][j-1]+match
                                    ,traceback,i,j);
                else 
                // if symbol at index i of 1st seq is mismatch with symbol at index j of 2nd seq
                // compare with 1) symbol at index i of 1st seq align with gap(insert to 2nd sequence)
                //              2) symbol at index j of 2nd seq align with gap(insert to 1st sequence)
                // and add traceback symbol at index [i][j]  
                    score[i][j] = max(score[i-1][j]+gap
                                    ,score[i][j-1]+gap
                                    ,score[i-1][j-1]+mismatch
                                    ,traceback,i,j);
            }
        }
        System.out.println("scoring table");
        printMatrix(score);
        traceback[0][0] = "end";
        //print traceback table
        System.out.println("tracback table");
        for(String[] row:traceback)
        {
            for(String col:row)
            {
                System.out.print(col+"\t");
            }
            System.out.println("");
        }
        System.out.println("");
        System.out.println("Best Score = "+score[m][n]);
        System.out.println("for Global Alignment between"+"\n"+seq1+"\nand\n"+seq2);
        System.out.println("scoring scheme: match = "+match+"\n\t  mismatch = "+mismatch
                            +"\n\t  gap ="+gap);
        System.out.println("-------------------------------------------------");
        //print path from traceback table
        String one = "";
        String two = "";
        printPath(seq1,seq2,traceback,m,n,one,two);
    }
    
    void printPath(StringBuilder seq1, StringBuilder seq2
                        , String[][] matrix, int row, int col
                        , String one, String two) 
    {
        if (row == 0 && col == 0 || "end".equals(matrix[row][col]) || matrix[row][col] == null) {
            System.out.println("possible alignment " + count);
            count++;
            String tmp1=new StringBuilder(one).reverse().toString();
            String tmp2=new StringBuilder(two).reverse().toString();
            System.out.println(tmp1);
            for (int vb = 0; vb < one.length(); vb++) {
                System.out.print("| ");
            }
            System.out.println("");
            System.out.println(tmp2);
            System.out.println("--------------------------------------");
        } else {
            for (int x = 0; x < matrix[row][col].length(); x++) {
                String t = String.valueOf(matrix[row][col].charAt(x));
                if("D".equals(t)) {
                    printPath(seq1, seq2,matrix,row - 1, col - 1, one+(seq1.charAt(col - 1)),two+(seq2.charAt(row - 1)));
                } else if ("L".equals(t)) {
                    printPath(seq1, seq2,matrix,row, col - 1, one+(seq1.charAt(col - 1)), two+"-");
                } else if ("U".equals(t)) {
                    printPath(seq1, seq2,matrix,row - 1, col, one+"-", two+(seq2.charAt(row - 1)));
                }
            }
        }
    }
    
    int max(int up, int left, int diag, String[][] matrix,int i, int j)
    {
        if(up == left && up == diag) 
        {
            matrix[i][j] = "UDL";
            return up;
        }
        if(up > left && up > diag)
        {
            matrix[i][j] = "U";
            return up;
        }
        if(left > up && left > diag)
        {
            matrix[i][j] = "L";
            return left;
        }
        if(diag > up && diag > left)
        {
            matrix[i][j] = "D";
            return diag;
        }
        if((up >= left && up > diag) || (left >= up && left > diag))
        {
            matrix[i][j] = "LU";
            return up;
        }
        if((diag >= left && diag > up) || (left >= diag && left > up))
        {
            matrix[i][j] = "DL";
            return diag;
        }
        else
        {
            matrix[i][j] = "DU";
            return up;
        }
    }
    
    Object[] read(FileReader seq, FileReader scoring) throws IOException
    {
        BufferedReader br = new BufferedReader(seq);
        StringBuilder read1 = new StringBuilder(br.readLine());
        StringBuilder read2 = new StringBuilder(br.readLine());
        br = new BufferedReader(scoring);
        int[] scr = new int[3]; int i = 0;
        String line;
        while( (line = br.readLine()) != null)
        {
            scr[i] = Integer.parseInt(line);
            i += 1;
        }
        
        Object[] obj = {read1, read2, scr};
        return obj;
    }
    
    void printMatrix(int[][] matrix)
    {
        for(int[] row:matrix)
        {
            for(int col:row)
            {
                System.out.print(col+"\t");
            }
            System.out.println("");
        }
        System.out.println("");
    }
    
    public static void main(String[] args) throws FileNotFoundException, IOException 
    {   
        // input file is text file
        // in scoring file, 1st line is match score 2nd is mismatch score and 3rd is gap penalty
        // sequence file has 2 line, each line is 1st and 2nd sequence
        FileReader scr = new FileReader("scoring.txt");
        FileReader seq = new FileReader("sequence.txt");
        
        PairwiseAlignment PA = new PairwiseAlignment();
        PA.align(seq, scr);
        
    }
}
