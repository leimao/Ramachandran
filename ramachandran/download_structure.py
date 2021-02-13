from typing import List
import os
import requests
import json
import asyncio
import aiohttp
import time
import tqdm

sem = asyncio.Semaphore(20)


def query_x_ray_entries(maximum_resolution: float = 1.5) -> List[str]:

    # RCSB search API
    # https://search.rcsb.org/redoc/index.html
    # https://search.rcsb.org/index.html#search-api

    pdb_ids = []

    url = "https://search.rcsb.org/rcsbsearch/v1/query"

    # query = {
    #     "query": {
    #         "type": "terminal",
    #         "service": "text",
    #     },
    #     "request_options": {
    #         "return_all_hits": True
    #     },
    #     "return_type": "entry",
    # }

    query = {
        "query": {
            "type":
            "group",
            "logical_operator":
            "and",
            "nodes": [{
                "type": "terminal",
                "service": "text"
            }, {
                "type":
                "group",
                "logical_operator":
                "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "exptl.method",
                            "operator": "exact_match",
                            "value": "X-RAY DIFFRACTION"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "operator": "less_or_equal",
                            "value": maximum_resolution,
                            "attribute": "rcsb_entry_info.resolution_combined"
                        }
                    },
                ]
            }]
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True
        }
    }

    response = requests.post(url, data=json.dumps(query))
    response.raise_for_status()

    response_status = response.status_code

    text_content = response.text
    json_data = json.loads(text_content)

    # print(response_status)

    # print(json_data["total_count"])
    # print(type(json_data["result_set"]))

    # print(type(json_data))

    # print(json_data["result_set"][0])

    for result in json_data["result_set"]:
        pdb_ids.append(result["identifier"])

    return pdb_ids


# async def async_download_pdbx(pdb_id: str, download_dir: str,
#                               maximum_num_connections: int = 20) -> None:

#     pdb_url = "https://files.rcsb.org/download/{}.cif".format(pdb_id)
#     pdb_filepath = os.path.join(download_dir, "{}.cif".format(pdb_id))

#     connector = aiohttp.TCPConnector(limit_per_host=maximum_num_connections)
#     async with aiohttp.ClientSession(connector=connector) as session:
#         async with session.get(pdb_url, timeout=20) as response:
#             if response.status == 200:
#                 content = await response.read()
#                 with open(pdb_filepath, "wb") as fhand:
#                     fhand.write(content)
#             else:
#                 print("Unable to download PDBx {}!".format(pdb_id))

async def async_download_pdbx(pdb_id: str, download_dir: str, session: aiohttp.ClientSession, sem: asyncio.Semaphore) -> None:
    

    pdb_url = "https://files.rcsb.org/download/{}.cif".format(pdb_id)
    pdb_filepath = os.path.join(download_dir, "{}.cif".format(pdb_id))

    # This timeout is useless.
    # timeout = aiohttp.ClientTimeout(total=30)
    # A workaround solution:
    # https://stackoverflow.com/a/64686124
    # https://github.com/aio-libs/aiohttp/issues/3203
    timeout = aiohttp.ClientTimeout(total=None,sock_connect=30,sock_read=30)

    async with sem:
        async with session.get(pdb_url, timeout=timeout) as response:
            if response.status == 200:
                content = await response.read()
                with open(pdb_filepath, "wb") as fhand:
                    fhand.write(content)
                # print("PDBx {} downloaded.".format(pdb_id))
            else:
                print("Unable to download PDBx {}!".format(pdb_id))


async def async_download_pdbxs(pdb_ids: str,
                               download_dir: str,
                               maximum_num_connections: int = 5) -> None:

    if not os.path.exists(download_dir):
        os.mkdir(download_dir)
    
    tasks = []
    # https://pawelmhm.github.io/asyncio/python/aiohttp/2016/04/22/asyncio-aiohttp.html
    sem = asyncio.Semaphore(5)

    connector = aiohttp.TCPConnector(limit_per_host=maximum_num_connections)
    async with aiohttp.ClientSession(connector=connector) as session:
        for pdb_id in pdb_ids:
            task = asyncio.create_task(
                async_download_pdbx(
                    pdb_id=pdb_id,
                    download_dir=download_dir,
                    session=session,
                    sem=sem,
                    )
                )
            tasks.append(task)
        
        responses = [
            await f for f in tqdm.tqdm(asyncio.as_completed(tasks),
                                    total=len(tasks))
        ]

        # try:
        #     await asyncio.gather(*tasks)
        # except Exception as e:
        #     print(repr(e))

        # No progress bar.
        # await asyncio.gather(*tasks)


def download_x_ray_pdbxs(download_dir: str,
                         maximum_resolution: float = 1.5) -> None:

    print("Querying PDB IDs of maximum resolution {}...".format(
        maximum_resolution))
    pdb_ids = query_x_ray_entries(maximum_resolution=maximum_resolution)

    print("Downloading {} PDBxs...".format(len(pdb_ids)))
    time_start = time.time()
    # Python 3.7+
    # asyncio.run(async_download_pdbs(pdb_ids=pdb_ids, download_dir=download_dir))
    # Python 3.5-3.6
    loop = asyncio.get_event_loop()
    loop.run_until_complete(
        async_download_pdbxs(pdb_ids=pdb_ids, download_dir=download_dir))
    time_end = time.time()
    time_elapsed = time_end - time_start
    print("Download Time Elapsed: {:02d}:{:02d}:{:02d}".format(
        int(time_elapsed // 3600), int(time_elapsed % 3600 // 60),
        int(time_elapsed % 60 // 1)))
