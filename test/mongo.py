from pymongo import MongoClient

client = MongoClient()
print(client)

db_names = client.list_database_names()
print(db_names)


