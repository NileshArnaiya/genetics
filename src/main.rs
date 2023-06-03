use std::env;
use std::fs::File;
use std::io::{BufReader, Read};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use std::io::{BufRead, BufWriter, Write, Result};
use rayon::prelude::*;
use rand::seq::SliceRandom; //
use std::error::Error;
use std::path::Path;
use std::io::SeekFrom;


struct SNP {
    chromosome: String,
    name: String,
    position: u32,
    cases_ref_count: u32,
    cases_alt_count: u32,
    controls_ref_count: u32,
    controls_alt_count: u32,
}

impl SNP {
    fn new(
        chromosome: String,
        name: String,
        position: u32,
        cases_ref_count: u32,
        cases_alt_count: u32,
        controls_ref_count: u32,
        controls_alt_count: u32,
    ) -> Self {
        SNP {
            chromosome,
            name,
            position,
            cases_ref_count,
            cases_alt_count,
            controls_ref_count,
            controls_alt_count,
        }
    }
}

fn read_bed_file(filename: &str) -> Vec<u8> {
    let file = File::open(filename).unwrap();
    let mut reader = BufReader::new(file);

    // Discard the first 3 bytes
    let mut buffer = [0; 3];
    reader.read_exact(&mut buffer).unwrap();

    let mut data = Vec::new();
    let chunk_size = 64 * 1024 * 1024; // 64 MB chunk size
    let mut chunk = vec![0; chunk_size];
    loop {
        let bytes_read = reader.read(&mut chunk).unwrap();
        if bytes_read == 0 {
            break;
        }
        data.extend_from_slice(&chunk[0..bytes_read]);
    }

    data
}






fn compute_p_value(observed_chi2: f64, permuted_chi2: Vec<f64>) -> f64 {
    // Compute the p-value by comparing the observed χ² test statistic
    // to the permuted χ² test statistics

    // Your implementation here
    let count = permuted_chi2.iter().filter(|&&x| x >= observed_chi2).count();
    let p_value = (count + 1) as f64 / (permuted_chi2.len() as f64 + 1.0);

    p_value
    //0.0 // Placeholder
}


fn calculate_chi_square(snp: &SNP) -> f64 {
    let cases_total_count = snp.cases_ref_count + snp.cases_alt_count;
    let controls_total_count = snp.controls_ref_count + snp.controls_alt_count;
    let total_count = cases_total_count + controls_total_count;

    if cases_total_count == 0 || controls_total_count == 0 {
        return 0.0; // Return a default value when the count is zero
    }

    let cases_aaf = snp.cases_alt_count as f64 / cases_total_count as f64;
    let controls_aaf = snp.controls_alt_count as f64 / controls_total_count as f64;
    let aaf = total_count as f64 / (cases_total_count + controls_total_count) as f64;

    let cases_ref_expected = (cases_total_count as f64 * (1.0 - aaf)).round();
    let cases_alt_expected = (cases_total_count as f64 * aaf).round();
    let controls_ref_expected = (controls_total_count as f64 * (1.0 - aaf)).round();
    let controls_alt_expected = (controls_total_count as f64 * aaf).round();

    let chi_square = ((snp.cases_ref_count as f64 - cases_ref_expected).powi(2) / cases_ref_expected)
        + ((snp.cases_alt_count as f64 - cases_alt_expected).powi(2) / cases_alt_expected)
        + ((snp.controls_ref_count as f64 - controls_ref_expected).powi(2) / controls_ref_expected)
        + ((snp.controls_alt_count as f64 - controls_alt_expected).powi(2) / controls_alt_expected);

    chi_square
}

fn calculate_observed_chi_square(snps: &[SNP]) -> Vec<f64> {
    let mut observed_chi_square_values: Vec<f64> = Vec::new();

    for snp in snps {
        let cases_ref_count = snp.cases_ref_count;
        let cases_alt_count = snp.cases_alt_count;
        let controls_ref_count = snp.controls_ref_count;
        let controls_alt_count = snp.controls_alt_count;

        let observed_matrix = vec![
            vec![cases_ref_count, cases_alt_count],
            vec![controls_ref_count, controls_alt_count],
        ];

        let mut chi_square = 0.0;

        for i in 0..2 {
            for j in 0..2 {
                let observed = observed_matrix[i][j] as f64;
                let expected = (cases_ref_count + controls_ref_count) as f64 * (cases_alt_count + controls_alt_count) as f64
                    / (cases_ref_count + cases_alt_count + controls_ref_count + controls_alt_count) as f64;
                chi_square += (observed - expected).powi(2) / expected;
            }
        }

        observed_chi_square_values.push(chi_square);
    }

    observed_chi_square_values
}

fn compute_alternate_allele_frequency(genotype_data: &[u8],
    num_individuals: usize, 
    num_snps: usize) -> Vec<f32> {
   let bits_per_individual = 2; // two bits per individual
   let bits_per_byte = 8;
   let individuals_per_byte = bits_per_byte / bits_per_individual;

   let mut result: Vec<f32> = vec![0.0; num_snps];

   let bytes_needed = (num_individuals + individuals_per_byte - 1) / individuals_per_byte;

   for i in 0..num_snps {
       let mut alternate_allele_count = 0;
       let mut missing_count:u32 = 0;

       for j in 0..num_individuals {
           let idx = i * bytes_needed + (j / 4);
           let byte = genotype_data[idx];
           let tmp: u8 = byte.wrapping_shr(((j * bits_per_individual) % bits_per_byte as usize).try_into().unwrap()) & 3;
               
           let shift = bits_per_byte - bits_per_individual * (j % 4) - bits_per_individual;
           // let allele = (byte >> shift) & 0b11;
           let allele:u8 = tmp.reverse_bits().rotate_left(2);
           if allele == 2 {
              
               missing_count += 1;
           } if allele == 0 {
               alternate_allele_count += 2;
           
           } if allele == 1 {
               alternate_allele_count += 1;
           } 
       }

       let total_count = num_individuals as u32 - missing_count;
       let freq = alternate_allele_count as f32 / (2 * total_count) as f32;
       result.push(freq);

   }

   result
}





fn count_lines_in_file(filename: &str) -> usize {
    
    let path = Path::new(filename);
    println!("Path {}", path.display());
    let file = File::open(&path).unwrap();
    let mut reader = BufReader::new(file);
    let mut buffer = [0; 4096];
    let mut count = 0;
    loop {
        let n = reader.read(&mut buffer).unwrap();
        if n == 0 {
            break;
        }
        count += buffer[..n].iter().filter(|&&c| c == b'\n').count();
    }

    count
}

fn get_fam_contents(fam_file: &str, num_lines: usize) -> Result<Vec<u8>> {
    let file = File::open(fam_file)?;
    let reader = BufReader::new(file);

    let mut labels = Vec::with_capacity(num_lines);

    for line in reader.lines().take(num_lines) {
        let line = line?;
        let fields: Vec<&str> = line.split_whitespace().collect();

        if fields.len() >= 6 {
            let label = fields[5].parse::<u8>().unwrap();
            labels.push(label);
        }
    }

    Ok(labels)
}

fn count_genotypes(v: &[u8], labels: &[usize], n: usize, m: usize) -> Result<Vec<Vec<f64>>> {
    let case: usize = 2;
    let byte_per_snp: usize = (n + 3) / 4;
    
    
    let counts = (0..m).into_par_iter().map(|i| {
        let mut case_aac: f64 = 0.0;
        let mut case_rac: f64 = 0.0;
        let mut control_aac: f64 = 0.0;
        let mut control_rac: f64 = 0.0;

        for j in 0..n {
            let index = i * byte_per_snp + j / 4;
            let offset = 2 * (j % 4);
            let byte = v.get(index).copied().unwrap_or(0u8).reverse_bits();
            let mask: u8 = 0b11000000 >> offset;
            let bits = (byte & mask) >> (6 - offset);
            if j < labels.len() {
            let label = labels[j];

            if label == case {
                if bits == 0b00 {
                    case_rac += 2.0;
                } else if bits == 0b11 {
                    case_aac += 2.0;
                } else if bits == 0b01 {
                    case_aac += 1.0;
                    case_rac += 1.0;
                }
            } else {
                if bits == 0b00 {
                    control_rac += 2.0;
                } else if bits == 0b11 {
                    control_aac += 2.0;
                } else if bits == 0b01 {
                    control_aac += 1.0;
                    control_rac += 1.0;
                }
            }
        }
        }

        vec![case_rac, case_aac, control_rac, control_aac]
    }).collect::<Vec<_>>();

    Ok(counts)
}



fn extract_labels_from_file(file_path: &str) -> std::io::Result<Vec<usize>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let labels: Vec<usize> = reader
        .lines()
        .map(|line| {
            let unwrapped_line = line.unwrap();
            let fields: Vec<&str> = unwrapped_line.split_whitespace().collect();
            fields[5].parse::<usize>().unwrap()
        })
        .collect();

    Ok(labels)
}


fn get_expected(aaf: &[f64], counts: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let mut adjusted_counts: Vec<Vec<f64>> = Vec::new();
    let len_diff = aaf.len() as isize - counts.len() as isize;

    if len_diff > 0 {
        // Pad the counts vector with zeros
        adjusted_counts.extend_from_slice(counts);
        for _ in 0..len_diff {
            adjusted_counts.push(vec![0.0; 2]);
        }
    } else if len_diff < 0 {
        // Truncate the counts vector
        adjusted_counts.extend_from_slice(&counts[..aaf.len()]);
    } else {
        adjusted_counts.extend_from_slice(counts);
    }

    let results: Vec<Vec<f64>> = aaf
        .iter()
        .zip(adjusted_counts)
        .map(|(&a, counts)| {
            let raf = 1.0 - a;
            vec![
                a * counts[0] + raf * counts[1],
            ]
        })
        .collect();

    Ok(results)
}

//
    // E = [ 
    //     [0.9999999999998, 9.0, 0.9999999998, 9.0],
    //     [4.0,6.0,4.0,6.0],
    //     [1.0,3.0, 0.0, 0.0],
    //     [3.428571428571429,4.571428571428571, 2.5714285,3.428571428571],
    //     [1.5,4.5,1.5,4.5]
    // ]

//

fn get_case_control_counts(real_data: &[u8], fam_labels: &[usize], n: usize, m: usize) -> Vec<Vec<f64>> {
    let mut counts = vec![vec![0.0; 4]; m];

    for (i, &label) in fam_labels.iter().enumerate() {
        let mut case_aac = 0.0;
        let mut case_rac = 0.0;
        let mut control_aac = 0.0;
        let mut control_rac = 0.0;

        for j in 0..n {
            let idx = i * n + j;
            let freq = real_data.get(idx).copied().unwrap_or(0);

            if label == 1 {
                if freq == 0 {
                    control_rac += 2.0;
                } else if freq == 1 {
                    control_aac += 2.0;
                } else if freq == 3 {
                    control_aac += 1.0;
                    control_rac += 1.0;
                }
            } else if label == 2 {
                if freq == 0 {
                    case_rac += 2.0;
                } else if freq == 1 {
                    case_aac += 2.0;
                } else if freq == 3 {
                    case_aac += 1.0;
                    case_rac += 1.0;
                }
            }
        }

        counts[i][0] = case_rac;
        counts[i][1] = case_aac;
        counts[i][2] = control_rac;
        counts[i][3] = control_aac;
    }

    counts
}

fn chi_square_test(observed: &[Vec<f64>], expected: &[Vec<f64>]) -> Result<Vec<f64>> {
    if observed.len() != expected.len() {
        println!("Gabdad hai");
    }

    let mut chi_square_values: Vec<f64> = Vec::new();

    for i in 0..observed.len() {
        let mut chi_square: f64 = 0.0;

        for j in 0..observed[i].len() {
            let o = observed[i][j];
            let e = expected[i][j];

            if e != 0.0 {
                chi_square += ((o - e).powi(2)) / e;
            }
        }

        chi_square_values.push(chi_square);
    }

    Ok(chi_square_values)
}

fn calculate_p_values(chi_square_values: &[f64]) -> Vec<f64> {
    let degrees_of_freedom = 1; // Degrees of freedom for the chi-square distribution
    
    let p_values: Vec<f64> = chi_square_values
        .iter()
        .map(|&chi_square| {
            let mut p_value = 0.0;
            if chi_square > 0.0 {
                let mut df = degrees_of_freedom as f64;
                let mut term = chi_square / 2.0;
                let mut sum = term.exp();
                
                for i in 1..(degrees_of_freedom / 2) {
                    term *= chi_square / (2.0 * df);
                    sum += term.exp();
                    df -= 1.0;
                }
                
                p_value = 1.0 - sum;
            }
            
            p_value
        })
        .collect();
    
    p_values
}


fn shuffle_data<T>(data: &mut [T]) {
    let mut rng = rand::thread_rng();
    for i in (1..data.len()).rev() {
        let j = rng.gen_range(0..=i);
        data.swap(i, j);
    }
}

fn print_type<T>(_: &T) {
    println!("{}", std::any::type_name::<T>());
}


fn main() {
    // Get command-line arguments
    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        panic!("Usage: permute <plink_base_name> <n>");
    }
    let mut plink_base_name = &args[1];
    let mut no_permutations = args[2].parse::<u32>().unwrap();

    println!("PLINK base name: {}.fam", plink_base_name);
    println!("Number of permutations: {}", no_permutations);

    // Read BIM file to get SNP information
    let mut bim_filename = format!("{}.bim", plink_base_name);
    let mut bim_file = match File::open(bim_filename.clone()) {
        Ok(file) => file,
        Err(_) => {
            println!("Failed to open BIM file");
            return;
        }
    };
    let bim_reader = BufReader::new(bim_file);
    let n = count_lines_in_file(&bim_filename.clone());
    println!("N(Bim) is {} ", n);
    let mut fam_file = format!("{}.fam", plink_base_name);
    let mut fam_contents = std::fs::read_to_string(&fam_file).expect("Error reading fam file");
    
    let m = count_lines_in_file(&fam_file);
    println!("M(Fam) is {} ", m);
    println!("Fam content: {}", fam_contents);

    let case_control: Vec<u8> = fam_contents.lines()
    .map(|line| {
        let fields: Vec<&str> = line.split_whitespace().collect();
        fields[5].parse::<u8>().expect("Error parsing case/control status")
    })
    .collect();
    println!("Case Control : {:?}", case_control);

    let real_data = read_bed_file(&(plink_base_name.to_owned() + ".bed"));
    println!("Real Bed : {:?}", real_data);
  
    let aaf = compute_alternate_allele_frequency(&real_data,m,n);
    println!("AAF are {:? }", aaf);

    let fam_labels = get_fam_contents(&fam_file, n).unwrap();
    println!("FAM are {:? }", fam_labels);

    //check
    let labels = match extract_labels_from_file(&fam_file) {
        Ok(labels) => labels,
        Err(err) => {
            panic!("Error loading labels: {}", err);
        }
    };

    // let fam_labels2 = read_fam_labels(fam_file_path)?;
    // println!("Fam Labels {} ", fam_labels2);
    // let genotype_data = read_bed_data(bed_file_path)?;
    // println!("Bed Data {} ", genotype_data);
    // let num_individuals = fam_labels2.len();
    // println!("Individuaks Data {} ", num_individuals);
    // let case_control: Vec<usize> = fam_labels2.iter().map(|&label| if label == 1 { 1 } else { 2 }).collect();

    // println!("Case contrrrrrr {} ", case_control);

    
    // read bim file --> get vector leng 
    // read fam file --> get vector leng (Case Controls)
    // read bed file 

        // we don't need first 3 bytes of bed file 
        // aaf add some paralleism 
        // 

        // for expected - loop thru each snp.. 2 * no.of case labels = case_count
        // labels.iter().filter(pass closure).count()  --> proprotion of total cases/count




    let fam_labels_usize: Vec<usize> = fam_labels.iter().map(|&label| label as usize).collect();
    println!("Fam Labels: {:?}", fam_labels_usize);
// let totals = bro_totals(&real_data, fam_labels_usize.as_slice(), n, m);

    let totals = get_case_control_counts(&real_data, fam_labels_usize.as_slice(), n, m);
    println!("Total Cases: {:?}", totals);

    let aaf_64 : Vec<f64> = aaf.iter().map(|&value| value as f64).collect();
    // let e = get_expected(&aaf_64, &totals);
    // println!("Lets see expected {:?}", e);
    
    let o: Vec<_> = count_genotypes(&real_data, &fam_labels_usize, n, m).unwrap();
    println!("Onserved is {:?}", o);

    let aaf: Vec<f64> = vec![0.0, 0.44444445, 0.5, 0.42857143, 0.25];
    
    let expected = get_expected(&aaf_64, &totals);

    println!("Expected: {:?}", expected);
    let observed: &[Vec<f64>] = &[
        vec![2.0, 8.0, 0.0, 10.0],
        vec![0.0, 10.0, 8.0, 2.0],
        vec![1.0, 3.0, 0.0, 0.0],
        vec![3.0, 5.0, 3.0, 3.0],
        vec![2.0, 4.0, 1.0, 5.0],
    ];

    let expected: &[Vec<f64>] = &[
        vec![0.9999999999998, 9.0, 0.9999999998, 9.0],
        vec![4.0, 6.0, 4.0, 6.0],
        vec![1.0, 3.0, 0.0, 0.0],
        vec![3.428571428571429, 4.571428571428571, 2.5714285, 3.428571428571],
        vec![1.5, 4.5, 1.5, 4.5],
    ];

    
    let chi_square_values = chi_square_test(&o, &totals).unwrap();
    println!("Results1 {:?}", chi_square_values);

    let result = chi_square_test(observed, expected).unwrap();
    println!("Results2 {:?}", result);

    let p1 = calculate_p_values(&chi_square_values);
    println!("P1 values {:?}", p1);


    let p = calculate_p_values(&result);
    println!("P values {:?}", p);

   let mut counter = 0;

    for _ in 0..no_permutations {
        let mut permuted_case_control = totals.to_vec();
        shuffle_data(&mut permuted_case_control);
        println!("Permuted {:?} ",permuted_case_control);

        let permuted_observed_counts: Vec<_> = count_genotypes(&real_data, &fam_labels_usize, n, m).unwrap();
        println!("Onserved is {:?}", o);

        // let permuted_expected_counts = get_expected(&aaf_64, &permuted_case_control);
        let permuted_chi_square = chi_square_test(&permuted_observed_counts, &permuted_case_control);

        // println!("Vector: {:?}", vec);
        if let Ok(permuted_vec) = permuted_chi_square.as_ref().map_err(|_| ()) {
            // Compare the values
            for (real_val, permuted_val) in chi_square_values.iter().zip(permuted_vec.iter()) {
                if *real_val > *permuted_val {
                    // Perform the desired action when real_chi_square > permuted_chi_square
                    // println!("this chi square is greater: {} > {}", real_val, permuted_val);
                    counter += 1;
                }
            }
        } else {
            println!("Fail");
        }

        // print_type(&chi_square_values);
        // print_type(&permuted_chi_square);
      
    }

    let p_value = (counter as f64) / (n as f64);
    println!("P-value: {}", p_value);

    let mut file = File::create("output.wassoc").expect("Failed to create file");

    writeln!(file, "{:<10} {}", "CHISQ", "P").expect("Failed to write to file");
 
    // Writing data
    for (chi_square, p_value) in chi_square_values.iter().zip(p.iter()) {
        writeln!(file, "{:<10} {}", chi_square, p_value).expect("Failed to write to file");
    }
 
    println!("Data written to 'output.wassoc' file.");

    
}