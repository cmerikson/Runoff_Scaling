import csv
import requests
import time

# Function to check if a dam is mentioned in the search results
def check_dam(items):
    for item in items:
        if "dam" in item["snippet"].lower():
            return "Yes"
    return "No"

# Read the input CSV file
input_file = r"C:\Users\cmeri\OneDrive - Dartmouth College\Research\Runoff_Scaling\Data\Dams\ReferenceStations_list.csv"
output_file = r"C:\Users\cmeri\OneDrive - Dartmouth College\Research\Runoff_Scaling\Data\Dams\DammedStations_list.csv"

# API key for Google Custom Search JSON API
api_key = "AIzaSyDdKmXTuQ12F5_ZFynYU66v-7b4JcCXjPk"
# Search engine ID for your custom search engine
search_engine_id = "c7f83bb5e108041a9"

# URL for the Google Custom Search JSON API
api_url = "https://www.googleapis.com/customsearch/v1"

with open(input_file, "r") as file:
    reader = csv.DictReader(file)
    headers = reader.fieldnames

    # Add new column header
    headers.append("Dam")

    # Create a list to hold the updated rows
    updated_rows = [headers]

    # Process each row in the CSV file
    for row in reader:
        river_name = row["river"]
        search_query = f'"{river_name}" and "dam"'

        # Parameters for the API request
        params = {
            "key": api_key,
            "cx": search_engine_id,
            "q": search_query,
            "num": 5
        }

        # Send the API request
        response = requests.get(api_url, params=params)
        data = response.json()

        # Check if a dam is mentioned in the search results
        dam_info = check_dam(data.get("items", []))

        # Append the updated row to the list
        row["Dam"] = dam_info
        updated_rows.append(list(row.values()))

        print(f"Search for '{river_name}' completed.")

        # Delay before the next search to avoid rate limits
        time.sleep(2)

# Write the updated rows to the output CSV file
with open(output_file, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerows(updated_rows)

print("CSV file updated successfully")
