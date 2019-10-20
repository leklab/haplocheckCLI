import java.util.Collection;
import java.io.File;
import java.util.HashMap;
import java.util.ArrayList;
import java.io.FileWriter;
import java.io.IOException;

import com.google.common.math.Quantiles;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;

import phylotree.Phylotree;
import phylotree.PhylotreeManager;

import contamination.VariantSplitter;
import importer.VcfImporter;
import vcf.Sample;

import contamination.HaplogroupClassifier;
import core.SampleFile;

import contamination.ContaminationDetection;
import contamination.ContaminationDetection.Status;
import contamination.objects.ContaminationObject;


public class haplocheck_contam{

	Phylotree phylotree;
	Collection<File> vcf_list;

	public haplocheck_contam(String directoryName){
		phylotree = PhylotreeManager.getInstance().getPhylotree("phylotree17.xml", "weights17.txt");
		/*	public static Collection<File> getVcfFiles(String directoryName) {
		File directory = new File(directoryName);
		return FileUtils.listFiles(directory, new WildcardFileFilter(Arrays.asList("*.vcf.gz", "*.vcf")), null);
		}
		*/
		vcf_list = Utils.getVcfFiles(directoryName);
	}

	public int run(){

		if (vcf_list.size() > 1){
        	System.out.println("Only one VCF file");
			return 1;
		}


		try{
			File file = vcf_list.iterator().next();


			VariantSplitter splitter = new VariantSplitter();
			VcfImporter reader = new VcfImporter();



			HashMap<String, Sample> mutationServerSamples = reader.load(file, false);

			ArrayList<String> profiles = splitter.split(mutationServerSamples);


			HaplogroupClassifier classifier = new HaplogroupClassifier();

			SampleFile haplogrepSamples = classifier.calculateHaplogrops(phylotree, profiles);


			ContaminationDetection contamination = new ContaminationDetection();
			ArrayList<ContaminationObject> result = contamination.detect(mutationServerSamples,
					haplogrepSamples.getTestSamples());


			contamination.writeReport("output", result);
			contamination.writeReportAsJson("output_json", result);
			writeSummary("output_summary", result);


		}catch(Exception e){
			return 1;
		}


		return 0;
	}


	private void writeSummary(String outSummary, ArrayList<ContaminationObject> contaminationList) throws IOException {
		int countYes = 0;
		int countNo = 0;
		ArrayList<Integer> distanceList = new ArrayList<Integer>();

		for (ContaminationObject cont : contaminationList) {

			if (cont.getStatus() == Status.YES) {
				countYes++;
				distanceList.add(cont.getDistance());
			} else if (cont.getStatus() == Status.NO) {
				countNo++;
			}
		}

		JsonObject result = new JsonObject();
		result.add("Yes", new JsonPrimitive(countYes));
		result.add("No", new JsonPrimitive(countNo));
		result.add("Distance", new JsonPrimitive(0.0));
		result.add("25Percentile", new JsonPrimitive(0.0));
		result.add("75Percentile", new JsonPrimitive(0.0));

		if (distanceList.size() > 0) {
			double distanceMedian = com.google.common.math.Quantiles.median().compute(distanceList);
			double percentile25 = Quantiles.percentiles().index(25).compute(distanceList);
			double percentile75 = Quantiles.percentiles().index(75).compute(distanceList);
			result.add("Distance", new JsonPrimitive(distanceMedian));
			result.add("25Percentile", new JsonPrimitive(percentile25));
			result.add("75Percentile", new JsonPrimitive(percentile75));
		}

		FileWriter wr = new FileWriter(outSummary);
		wr.write(result.toString());
		wr.close();
	}



    public static void main(String[] args) {
        if(args.length != 1){
          System.out.println("Usage: java -jar mtServerCLI.jar bam_file");
          System.exit(1);
        }


        haplocheck_contam pileup = new haplocheck_contam(args[0]);
        //pileup.run();

    }



}