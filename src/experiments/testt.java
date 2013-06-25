package experiments;

public class testt {
	public static void main(String[] args) {
		
		int i = 0, j = 0, k =0;
		int N = 8;
		for (; i<N; i++){
			for(j=0;j<i+1;j++){
				System.out.print("*");
				for(k=0;k<i;k++){
					System.out.print(".");
				}
			}
			System.out.println("");
		}
	}
}
