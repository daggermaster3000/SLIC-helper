from tabulate import tabulate
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

def primers_to_df(inserts):
    headers = ["Insert", "", "Primer Type", "Sequence 5'->3'", "Tm", "GC Content", "Start Position"]
    table_data = []

    for insert, data in inserts.items():
        fw_primer = data.get("FW_primer", ["Not found"])
        rv_primer = data.get("RV_primer", ["Not found"])
        table_data.append([insert, "", "Forward Primer", *fw_primer])
        table_data.append(["", "", "Reverse Primer", *rv_primer])
    df = pd.DataFrame(table_data, columns=headers)
   

    return df

def read_seqs(file_path):

    # Dictionary to store the sequences
    sequences = {}

    # Temporary variables for keys and values
    current_key = None
    current_sequence = Seq("")

    # Open and read the file line by line
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                pass
            else:
                if line.startswith('>'):  # New sequence name detected
                    # Save the previous sequence if it exists
                    if current_key is not None:
                        
                        sequences[current_key] = {"Sequence":Seq("").join(current_sequence),"FW_primer":None,"RV_primer":None}
                    
                    # Reset the current sequence and set the new key
                    current_key = line[1:].strip()  # Remove '>' and trim any whitespace
                    current_sequence = []
                else:
                    # Add this line to the current sequence, removing unwanted characters like dots
                    current_sequence.append(line.strip())

    # Save the last sequence
    if current_key is not None:
        sequences[current_key] = {"Sequence":Seq("").join(current_sequence),"FW_primer":None,"RV_primer":None}

    # Print the result
    for key, value in sequences.items():
        if '*' in key:
            sequences[key]["Sequence"] = sequences[key]["Sequence"].reverse_complement()
          
    return sequences

def find_primer(sequence, max_primer_start_pos, window_size, max_primer_length, m_t, min_gc):
    """
    Find primers by widening a window then sliding it 
    """
    for i in range(max_primer_start_pos):
        
        for pl in range(window_size, max_primer_length):
            primer = sequence[i:i+pl]
            primer_tm = mt.Tm_GC(primer)
            gc_content = gc_fraction(primer)
            
            if primer_tm > m_t and gc_content > min_gc :
                
                if i != 0:
                    primer = sequence[0:i]+"-"+primer
                  
                    return  [primer, primer_tm, gc_content, i] 
                
                return  [primer, primer_tm, gc_content, i] # Exit the function when primer found
            

    return ["Not found","Not found","Not found"]


def get_primers(inserts,window_size=10,max_primer_length=30,max_primer_start_pos=10,min_gc=0.5,m_t=55):

    for insert,data in inserts.items():
        data["FW_primer"] = None
        data["RV_primer"] = None

        # fw primer
        gc_adj = 0
        m_t_adj = 0
        data["FW_primer"] = find_primer(data['Sequence'],max_primer_start_pos,window_size,max_primer_length,m_t,min_gc)
        
        while "Not found" in data["FW_primer"]:
            data["FW_primer"] = find_primer(data['Sequence'],max_primer_start_pos,window_size,max_primer_length,m_t-m_t_adj,min_gc-gc_adj)
            
            
            if min_gc-gc_adj > 0.1:
                gc_adj += 0.1
                
            if m_t-m_t_adj > 45:
                m_t_adj += 0.1
                
            elif data["FW_primer"][0] == "Not found":
                data["FW_primer"] = ["Not Found."]

        # rv primer
        gc_adj = 0
        m_t_adj = 0
        data["RV_primer"] = find_primer(data['Sequence'].reverse_complement(),max_primer_start_pos,window_size,max_primer_length,m_t,min_gc)
        while "Not found" in data["RV_primer"]:
            data["RV_primer"] = find_primer(data['Sequence'].reverse_complement(),max_primer_start_pos,window_size,max_primer_length,m_t-m_t_adj,min_gc-gc_adj)
                
            if min_gc-gc_adj > 0.1:
                gc_adj += 0.1
            if m_t-m_t_adj > 45:
                m_t_adj += 0.1
            elif data["RV_primer"][0] == "Not found":
                data["RV_primer"] = ["Not Found."]
    return inserts

def add_homology(inserts,homology_seq_length=10):
    insert_names = list(inserts.keys())
    for i,name in enumerate(insert_names):
        
        if i == len(insert_names)-1:
            hom_seq = inserts[insert_names[0]]["Sequence"].complement()[:homology_seq_length]
            inserts[insert_names[i]]["RV_primer"][0] = f'{hom_seq[::-1]}-{inserts[insert_names[i]]["RV_primer"][0]}'
        else:
            hom_seq = inserts[insert_names[i+1]]["Sequence"].complement()[:homology_seq_length]
            inserts[insert_names[i]]["RV_primer"][0] = f'{hom_seq[::-1]}-{inserts[insert_names[i]]["RV_primer"][0]}'

        if i == 0:
            hom_seq = inserts[insert_names[-1]]["Sequence"][-homology_seq_length:]
            inserts[insert_names[i]]["FW_primer"][0] = f'{hom_seq}-{inserts[insert_names[i]]["FW_primer"][0]}'
        else:
            hom_seq = inserts[insert_names[i-1]]["Sequence"][-homology_seq_length:]
            inserts[insert_names[i]]["FW_primer"][0] = f'{hom_seq }-{inserts[insert_names[i]]["FW_primer"][0]}'
            
    return inserts

def export(inserts,homology_seq_length=10):
    display_length = 50
    insert_names = list(inserts.keys())
    for id, name in enumerate(insert_names):

        # append primers to sequence for visualization purposes
        if id == 0:
            poi = insert_names[id]  # part of interest
            homology =  "".join(sorted([poi,insert_names[-1]])).replace("*","rev")
            inserts[poi]["FW_sequence_woverhang"]=f'<span class="{homology}">{inserts[poi]["FW_primer"][0].split("-")[0]}</span>{inserts[poi]["Sequence"][:display_length]}...{inserts[poi]["Sequence"][-display_length:]}{"-"*homology_seq_length}'
        
        else:
            poi = insert_names[id]  # part of interest
            homology =  "".join(sorted([poi,insert_names[id-1]])).replace("*","rev")
            inserts[poi]["FW_sequence_woverhang"]=f'<span class="{homology}">{inserts[poi]["FW_primer"][0].split("-")[0]}</span>{inserts[poi]["Sequence"][:display_length]}...{inserts[poi]["Sequence"][-display_length:]}{"-"*homology_seq_length}'
            
        if id == len(insert_names)-1:
            poi = insert_names[id]  # part of interest
            
            homology =  "".join(sorted([poi,insert_names[0]])).replace("*","rev")
            inserts[poi]["RV_sequence_woverhang"]=f'{"-"*homology_seq_length}{inserts[poi]["Sequence"].complement()[:display_length]}...{inserts[poi]["Sequence"].complement()[-display_length:]}<span class="{homology}">{inserts[poi]["RV_primer"][0].split("-")[0][::-1]}</span>'
        else:
            poi = insert_names[id]  # part of interest
            homology =  "".join(sorted([poi,insert_names[id+1]])).replace("*","rev")
            inserts[poi]["RV_sequence_woverhang"]=f'{"-"*homology_seq_length}{inserts[poi]["Sequence"].complement()[:display_length]}...{inserts[poi]["Sequence"].complement()[-display_length:]}<span class="{homology}">{inserts[poi]["RV_primer"][0].split("-")[0][::-1]}</span>'

# Create final construct
def create_final_construct(inserts):
    finalconstruct = Seq("")
    outname = ""
    for insert in inserts:
        outname += str(insert)
        finalconstruct += inserts[insert]["Sequence"]
    return finalconstruct

js = """
<script>
// Define a function to generate a random color
function getRandomColor() {
  const letters = '0123456789ABCDEF';
  let color = '#';
  for (let i = 0; i < 6; i++) {
    color += letters[Math.floor(Math.random() * 16)];
  }
  return color;
}

// Get all elements in the DOM
const allElements = document.querySelectorAll('*');

// Create a Set to store unique class names
const uniqueClasses = new Set();

// Iterate through all elements to collect unique class names
allElements.forEach(element => {
  const classes = element.classList;
  classes.forEach(className => {
    if (className != 'dataframe'){
    uniqueClasses.add(className);
    }
  });
});

// Assign a random color to all elements sharing each class
uniqueClasses.forEach(className => {
  const elementsWithClass = document.querySelectorAll('.' + className);
  const randomColor = getRandomColor();
  elementsWithClass.forEach(element => {
    element.style.backgroundColor = randomColor;
  });
});
</script>
<style>
p {
    font-family: Consolas, Monaco, 'Andale Mono', 'Ubuntu Mono', monospace;
    font-size: 12px;
    line-height: 1.5;
    padding: 10px;
    background-color: #f4f4f4;
    border: 1px solid #ddd;
    border-radius: 5px;
    overflow-x: auto; /* Enable horizontal scrolling if needed */
    white-space: pre-wrap; /* Preserve line breaks */
}
/* Reset default table styles */
table {
  border-collapse: collapse;
  width: 100%;
}

/* Table header styles */
thead {
  background-color: #f2f2f2 !important;
}

/* Table cell styles */
td, th {
  font-family: Consolas, Monaco, 'Andale Mono', 'Ubuntu Mono', monospace;
  border: 1px solid #dddddd;
  text-align: left;
  padding: 8px;
  font-size: 12px;
}

/* Alternate row color */


/* Hover effect */
tbody tr:hover {
  background-color: #f2f2f2 !important;
}

/* Modern table style */
.modern-table {
  border-radius: 8px;
  overflow: hidden;
  box-shadow: 0 0 20px rgba(0, 0, 0, 0.1);
}
h2{
font-family: Consolas, Monaco, 'Andale Mono', 'Ubuntu Mono', monospace;
}
h3{
font-family: Consolas, Monaco, 'Andale Mono', 'Ubuntu Mono', monospace;
font-size: 15px;
}



</style>
"""
Title="""
<div style="display:ruby;"><h2>SLIC helper 🧬</h2><h3 >v1.0</h3></div>
<h3>Notes:</h3>
<p>Below are the primers found by the algorithm in [homology-ext-binding] format.<br>The Tm is for the binding sequence only<br>Us ctrl-f to verify the results</p>
<h3>Primers:</h2>
"""

if __name__ == '__main__' :
    inserts = read_seqs('/Users/quillan/Documents/Lab/Thesis/Random stuff/primerstuff/QF_Pkd2l1_Inpp5e.txt')
    outpath = "output"
    
    # find primers
    inserts = get_primers(inserts)

    # append homology sequence overhang
    inserts = add_homology(inserts)

    # convert to df
    primers_df = primers_to_df(inserts)
    html_df = primers_df.to_html()

    # add export info to the inserts
    export(inserts)

    # create final construct
    finalconstruct = create_final_construct(inserts)
    with open(outpath+".html", 'w') as file:
            file.write(Title)
            file.write(html_df)
            for insert,data in inserts.items():
                    file.write(f'<h3>{insert}</h3>')
                    file.write(f'<p>{data["FW_sequence_woverhang"]}<br>{data["RV_sequence_woverhang"]}</p>')
                    #file.write(f'{data["RV_sequence_woverhang"]}')
            file.write(f'<h2>Final construct</h2>')
            file.write(f'<p>{str(finalconstruct)}<br>{str(finalconstruct.complement())}</p>') 
            #file.write(f'<p>{str(finalconstruct.complement())}</p>')     
            file.write(js)