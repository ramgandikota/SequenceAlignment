
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Scanner;

public class SequenceAlignment {
	
	String alphabetStr = null;
	int alphabetLen = 0;
	HashMap<Character,Integer> alphabetMap = new HashMap<>();
	StringBuilder sb = null;
	ArrayList<String> scoringList = new ArrayList<>();
	ArrayList<String> queryList = new ArrayList<>();
	HashMap<String, Integer> queryIdMap = new HashMap<>();
	HashMap<String, Integer> databaseIdMap = new HashMap<>();
	ArrayList<String> databaseList = new ArrayList<>();
	static int [][] scoringMatrix; 
	
	

public void readAll(String queryFilePath,String databaseFilePath, String alphabetPath, String scorMatPath ) throws IOException{
	
	Scanner scan = null;
	
	//Reading the Alphabet 
scan = new Scanner(new File(alphabetPath));
alphabetStr = scan.nextLine();
alphabetLen = alphabetStr.length();
for(int i=0;i<alphabetLen;i++){
	alphabetMap.put(alphabetStr.toLowerCase().charAt(i), i);
}


//Reading and creating Scoring Matrix 

scoringMatrix = new int[alphabetLen][alphabetLen];	 
scan = new Scanner(new File(scorMatPath));

for( int i=0; i<alphabetLen ; i++){
	for( int j=0; j<alphabetLen;j++){
		
		scoringMatrix[i][j] = scan.nextInt();	
	}
	
 }

// Query Input reading and creating query list

scan = new Scanner(new File(queryFilePath));
sb = new StringBuilder();
String line = null;
int  id = 0;
 line = scan.nextLine();
 id = Integer.parseInt(line.split("\\s")[0].substring(5,line.split("\\s")[0].length()));
 while(scan.hasNext()){	
	 line = scan.nextLine();
	 if(!line.startsWith(">hsa")){
	 
		 sb.append(line);
	 
 }else{
	 if( sb.length() > 0){
		 
		 queryList.add(sb.toString());
		 queryIdMap.put(sb.toString(),id);
		 sb = new StringBuilder();
		 id = Integer.parseInt(line.split("\\s")[0].substring(5,line.split("\\s")[0].length()));
			 
		 }
		 else{
			 line = scan.nextLine();
			 sb.append(line);
			 
		 }

	 }
		
 }
	 queryList.add(sb.toString());
	 queryIdMap.put(sb.toString(),id);
	 
// Database reading and creating database list		
sb = new StringBuilder();
scan = new Scanner(new File(databaseFilePath));
 line = scan.nextLine();
 id = Integer.parseInt(line.split("\\s")[0].substring(5,line.split("\\s")[0].length()));
 while(scan.hasNext()){	
	 line = scan.nextLine();
	 if(!line.startsWith(">hsa")){
	 
	 sb.append(line);
 }else{
	 if( sb.length() > 0){
			
		 databaseList.add(sb.toString());
		 databaseIdMap.put(sb.toString(),id);
		 sb = new StringBuilder();
		 id = Integer.parseInt(line.split("\\s")[0].substring(5,line.split("\\s")[0].length()));
	 }else{

		 line = scan.nextLine();
		 sb.append(line);
		}
	}
}
 databaseList.add(sb.toString());
 databaseIdMap.put(sb.toString(),id);
 scan.close();
	
}

private Results globalAlignment(String seq1, String seq2, int gap) {
	// TODO Auto-generated method stub
int id1 = queryIdMap.get(seq1), id2 = databaseIdMap.get(seq2);
int[][] costMatrix = new int[seq1.length()+1][seq2.length()+1];

//Initialize the score matrix
//the first row and column are for the gap
//Complexity: O(NxM)

for (int i =0; i< seq1.length()+1; i++)
{
	for (int j =0; j< seq2.length()+1; j++)
	{
		costMatrix[i][j] =0;
		if (i==0){
			costMatrix[i][j] = gap*j;
		}else if (j==0) {
			costMatrix[i][j] = gap*i;
		}
	}
}
             
int similarityCost=0;
int matchCost=0;
int seq1GapCost=0;
int seq2GapCost=0;
             
//Compute the minimum cost scores between all 
//possible pairs of prefixes
//Complexity: O(NxM)
for (int i =1; i< seq1.length()+1; i++){
	for (int j =1; j< seq2.length()+1; j++){
		
		//Case 1: The cost of mismatch between the two prefixes

similarityCost= scoringMatrix[alphabetMap.get(seq1.charAt(i-1))][alphabetMap.get(seq2.charAt(j-1))];	
matchCost = costMatrix[i-1][j-1] + similarityCost;

//Case 2: the cost of adding a gap on sequence 2
seq2GapCost = costMatrix[i-1][j] + gap;
//Case 3: the cost of adding a gap on sequence 1
		seq1GapCost = costMatrix[i][j-1] + gap;
		costMatrix[i][j] = Math.max(Math.max(matchCost,seq1GapCost),
                             seq2GapCost);
 }
}        

//BackTracking
StringBuilder alignedSeq1= new StringBuilder();
StringBuilder alignedSeq2= new StringBuilder();
             
int j = seq2.length();
int i = seq1.length();
             
while (i >0 || j > 0) {
 if (i>0 && j > 0){
	 similarityCost= scoringMatrix[alphabetMap.get(seq1.charAt(i-1))][alphabetMap.get(seq2.charAt(j-1))];
 }
          
if ( i > 0 && j >0 && costMatrix[i][j] == costMatrix[i-1][j-1] + similarityCost) { 
    
	alignedSeq1.append(seq1.charAt(i-1));
    alignedSeq2.append(seq2.charAt(j-1));
    i=i-1;
    j=j-1;
 }
 else if ( i > 0 && costMatrix[i][j] == costMatrix[i-1][j] + gap){
    alignedSeq2.append("-");
    alignedSeq1.append(seq1.charAt(i-1));
    i=i-1;
 }
else if ( j > 0 && costMatrix[i][j] == costMatrix[i][j-1] + gap){
    alignedSeq1.append("-");
    alignedSeq2.append(seq2.charAt(j-1));
    j=j-1;
}
} // end of while
	 
	 return new Results(0, 0, costMatrix[seq1.length()][seq2.length()], new String[]{ alignedSeq1.reverse().toString(), alignedSeq2.reverse().toString() }, id1, id2);
	 
	 }


private Results localAlignment(String seq1, String seq2, int gap) {
	// TODO Auto-generated method stub
int id1 = queryIdMap.get(seq1), id2 = databaseIdMap.get(seq2);
int[][] costMatrix = new int[seq1.length()+1][seq2.length()+1];
int maxX=0, maxY=0, maxCost = Integer.MIN_VALUE;


int similarityCost=0;
int matchCost=0;
int seq1GapCost=0;
int seq2GapCost=0;
             
//Compute the minimum cost scores between all 
//possible pairs of prefixes
//Complexity: O(NxM)
for (int i =1; i<seq1.length()+1; i++){
	for (int j =1; j<seq2.length()+1; j++){
		
		//Case 1: The cost of mismatch between the two prefixes
similarityCost= scoringMatrix[alphabetMap.get(seq1.charAt(i-1))][alphabetMap.get(seq2.charAt(j-1))];	
matchCost = costMatrix[i-1][j-1] + similarityCost;

//Case 2: the cost of adding a gap on sequence 2
seq2GapCost = costMatrix[i-1][j] + gap;
//Case 3: the cost of adding a gap on sequence 1
		seq1GapCost = costMatrix[i][j-1] + gap;
		
		costMatrix[i][j]=Math.max(Math.max(0, seq2GapCost),Math.max(matchCost,seq1GapCost));
		if(costMatrix[i][j] > maxCost){
			maxCost = costMatrix[i][j];
			maxX= i;
			maxY =j;
		}
 }
}        

//BackTracking
StringBuilder alignedSeq1= new StringBuilder();
StringBuilder alignedSeq2= new StringBuilder();
             
int j = maxY;
int i = maxX;
             
while (i>0 && j>0 && costMatrix[i][j] != 0) {
	
 if (i>0 && j > 0){
	similarityCost= scoringMatrix[alphabetMap.get(seq1.charAt(i-1))][alphabetMap.get(seq2.charAt(j-1))];
	 //similarityCost = (seq2.charAt(j-1) == seq1.charAt(i-1))?5:-4;
 }
          
if (i > 0 && j >0 && costMatrix[i][j] == costMatrix[i-1][j-1] + similarityCost) { 
    
	alignedSeq1.append(seq1.charAt(i-1));
    alignedSeq2.append(seq2.charAt(j-1));
    i=i-1;
    j=j-1;
 }
 else if ( i> 0 && costMatrix[i][j] == costMatrix[i-1][j] + gap){
    alignedSeq2.append("-");
    alignedSeq1.append(seq1.charAt(i-1));
    i=i-1;
 }
else if ( j>0 && costMatrix[i][j] == costMatrix[i][j-1] + gap){
    alignedSeq1.append("-");
    alignedSeq2.append(seq2.charAt(j-1));
    j=j-1;
}
} // end of while
	 
	
	 return new Results(i, j, costMatrix[maxX][maxY], new String[]{ alignedSeq1.reverse().toString(), alignedSeq2.reverse().toString() }, id1, id2);
}

private Results doveTailAlignment(String seq1, String seq2, int gap) {
	// TODO Auto-generated method stub
int id1 = queryIdMap.get(seq1), id2 = databaseIdMap.get(seq2);
int[][] costMatrix = new int[seq1.length()+1][seq2.length()+1];
int row_maxX=0, row_maxY=0, col_maxX=0, col_maxY=0,maxCostRow = Integer.MIN_VALUE, maxCostColumn = Integer.MIN_VALUE;


//Initialize the score matrix
//the first row and column are for the gap
//Complexity: O(NxM)
for (int i =0; i< seq1.length()+1; i++)
{
	for (int j =0; j< seq2.length()+1; j++)
	{
		if (i==0){
			costMatrix[i][j] = 0;
		}else if (j==0) {
			costMatrix[i][j] = 0;
		}
	}
}
             
int similarityCost=0;
int matchCost=0;
int seq1GapCost=0;
int seq2GapCost=0;
             
//Compute the minimum cost scores between all 
//possible pairs of prefixes
//Complexity: O(NxM)
for (int i =1; i< seq1.length()+1; i++){
	for (int j =1; j< seq2.length()+1; j++){
		
		//Case 1: The cost of mismatch between the two prefixes
similarityCost= scoringMatrix[alphabetMap.get(seq1.charAt(i-1))][alphabetMap.get(seq2.charAt(j-1))];	
matchCost = costMatrix[i-1][j-1] + similarityCost;

//Case 2: the cost of adding a gap on sequence 2
seq2GapCost = costMatrix[i-1][j] + gap;
//Case 3: the cost of adding a gap on sequence 1
		seq1GapCost = costMatrix[i][j-1] + gap;
		costMatrix[i][j] = Math.max(Math.max(matchCost,seq1GapCost),seq2GapCost);
		
		if(i == seq1.length()){
			if(costMatrix[i][j] > maxCostRow){
				maxCostRow = costMatrix[i][j];
				row_maxX = i;
				row_maxY = j;
			}
		}else if(j == seq2.length()){
			if(costMatrix[i][j] > maxCostColumn){
				maxCostColumn = costMatrix[i][j];
				col_maxX = i;
				col_maxY = j;
			}
		}
		
 }
}        

//BackTracking
StringBuilder alignedSeq1= new StringBuilder();
StringBuilder alignedSeq2= new StringBuilder();
int i =0, j=0, startPos1 = 0, startPos2 =0;

  if(maxCostRow == Math.max(maxCostRow, maxCostColumn)){
	  i=row_maxX;
	  j=row_maxY;
	  
  }else{
	  i=col_maxX;
	  j=col_maxY;
	  
  }
  startPos1 = j;
  startPos2 = i;
while (i>0 && j>0 && costMatrix[i][j] != 0) {
	
	
 if (i>0 && j > 0){
	similarityCost= scoringMatrix[alphabetMap.get(seq1.charAt(i-1))][alphabetMap.get(seq2.charAt(j-1))];
 }
          
if ( i>0 && j>0 && costMatrix[i][j] == costMatrix[i-1][j-1] + similarityCost) { 
    
	alignedSeq1.append(seq1.charAt(i-1));
    alignedSeq2.append(seq2.charAt(j-1));
    i=i-1;
    j=j-1;
 }
 else if ( i>0 && costMatrix[i][j] == costMatrix[i-1][j] + gap){
    alignedSeq2.append("-");
    alignedSeq1.append(seq1.charAt(i-1));
    i=i-1;
 }
else if ( j>0 && costMatrix[i][j] == costMatrix[i][j-1] + gap){
    alignedSeq1.append("-");
    alignedSeq2.append(seq2.charAt(j-1));
    j=j-1;
}
} // end of while

	 return new Results(i, j, costMatrix[startPos2][startPos1], new String[]{ alignedSeq1.reverse().toString(), alignedSeq2.reverse().toString() }, id1, id2);
}



public static void main(String[] args) throws IOException {

	SequenceAlignment align = new SequenceAlignment(); 
	Results result = null;
	int k = 0;
	ArrayList<Results> resultList= new ArrayList<>();
	int numberOfResults = Integer.parseInt(args[5]);
	
	PriorityQueue<Results> minheap=new PriorityQueue<Results>(numberOfResults,new Comparator<Results>() {

		@Override
		public int compare(Results o1, Results o2) {
			// TODO Auto-generated method stub
		int res = o1.score-o2.score;
		return res;
	}  
});

//Reading all the files from the command line arguments
align.readAll(args[1],args[2],args[3],args[4]);

//Switch case to choose between Global, Local or Dovetail alignment
 int option = Integer.parseInt(args[0]);
 switch(option){
 case 1:
	 for(int i=0;i<align.queryList.size();i++){
		 for(int j=0;j<align.databaseList.size();j++){
			 
			 result = align.globalAlignment(align.queryList.get(i),align.databaseList.get(j),Integer.parseInt(args[6]));
			 minheap.add(result);
		 }
	 }
	 k = minheap.size();
	 for( int i=0; i<(k-numberOfResults);i++){				
		 minheap.remove();
	 }
	 break;
 case 2:
	 for(int i=0;i<align.queryList.size();i++){
		 for(int j=0;j<align.databaseList.size();j++){
			 result = align.localAlignment(align.queryList.get(i),align.databaseList.get(j),Integer.parseInt(args[6]));
			 minheap.add(result);
		 }
	 }
	 k = minheap.size();
	 for( int i=0; i<(k-numberOfResults);i++){
		 minheap.remove();
	 }
	 break;
 case 3:
	 for(int i=0;i<align.queryList.size();i++){
		 for(int j=0;j<align.databaseList.size();j++){
			 
			result =  align.doveTailAlignment(align.queryList.get(i),align.databaseList.get(j),Integer.parseInt(args[6]));
			 minheap.add(result);
		 }
	 }
	 k = minheap.size();
	 for( int i=0; i<(k-numberOfResults);i++){
		 minheap.remove();
	 }	 
	 break;
default:
	System.out.println("Invalid option, choose between 1, 2 and 3");
 }
 
//Removing the max elements from Heap and adding the results to a list 
for( int i=0 ; i<numberOfResults;i++){
	resultList.add(minheap.remove());
}

//Printing the max elements
for( int i=numberOfResults-1;i>=0;i--){
	
	System.out.println("Score="+resultList.get(i).score+"\n");
System.out.println(resultList.get(i).id1+" "+resultList.get(i).startPos1+" "+resultList.get(i).alignedSeq[0]);
System.out.println(resultList.get(i).id2+" "+resultList.get(i).startPos2+" "+resultList.get(i).alignedSeq[1]);
	}
 }
}
