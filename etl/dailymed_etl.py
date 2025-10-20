import zipfile
import os
import xmltodict
import pandas as pd

def extract_and_process_dailymed_data():
    """
    Extracts the DailyMed data and processes the XML files.
    """
    if not os.path.exists("dailymed_data"):
        os.makedirs("dailymed_data")

    with zipfile.ZipFile("dm_spl_release_human_rx_part1.zip", 'r') as zip_ref:
        zip_ref.extractall("dailymed_data")

    print("Extracted DailyMed data.")

    extracted_data = []

    # The zip file contains a 'prescription' directory with more zip files.
    # We need to unzip those as well.
    for root, dirs, files in os.walk("dailymed_data/prescription"):
        for filename in files:
            if filename.endswith(".zip"):
                zip_path = os.path.join(root, filename)
                with zipfile.ZipFile(zip_path, 'r') as inner_zip_ref:
                    for inner_filename in inner_zip_ref.namelist():
                        if inner_filename.endswith(".xml"):
                            with inner_zip_ref.open(inner_filename) as xml_file:
                                try:
                                    xml_string = xml_file.read()
                                    sections = parse_spl_document(xml_string)
                                    if sections:
                                        extracted_data.append(sections)
                                except Exception as e:
                                    print(f"Error processing {inner_filename}: {e}")

    df = pd.DataFrame(extracted_data)
    df.to_csv("dailymed_processed_data.csv", index=False)
    print("Processed DailyMed data and saved to dailymed_processed_data.csv")

def parse_spl_document(xml_string):
    """
    Parses the SPL XML and extracts the "ADVERSE REACTIONS", "WARNINGS", and "DOSAGE AND ADMINISTRATION" sections.
    """
    data_dict = xmltodict.parse(xml_string)

    sections = {}

    try:
        components = data_dict['document']['component']['structuredBody']['component']
        if not isinstance(components, list):
            components = [components]

        for component in components:
            section = component.get('section', {})
            if not section:
                continue

            title = section.get('title')
            if isinstance(title, dict):
                title = title.get('#text')

            if title and title.strip().upper() in ["ADVERSE REACTIONS", "WARNINGS AND PRECAUTIONS", "DOSAGE AND ADMINISTRATION"]:
                text_content = []
                text_element = section.get('text', {})

                if isinstance(text_element, dict):
                    paragraphs = text_element.get('paragraph', [])
                    if not isinstance(paragraphs, list):
                        paragraphs = [paragraphs]
                    for p in paragraphs:
                        if p: # Check if paragraph is not None
                            if isinstance(p, dict):
                                text_content.append(p.get('#text', ''))
                            else:
                                text_content.append(str(p))


                sections[title.strip()] = " ".join(text_content)

    except (KeyError, TypeError):
        pass # Handle cases where the structure is different

    return sections

if __name__ == "__main__":
    extract_and_process_dailymed_data()
