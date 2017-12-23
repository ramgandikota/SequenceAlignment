public class Results {
	
	String [] alignedSeq = null;
	int id1, id2, score,startPos1, startPos2;
	
	public Results(int startPos1, int startPos2, int score, String[] alignedSeq, int id1, int id2){
		this.startPos1 = startPos1;
		this.startPos2 = startPos2;
		this.score = score;
		this.alignedSeq = alignedSeq;
		this.id1 = id1;
		this.id2 = id2;
	}
	
}
