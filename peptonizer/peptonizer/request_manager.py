import requests


class RequestManager:
    @staticmethod
    def perform_post_request(url, payload):
        req = requests.Request('POST', url=url, json=payload)
        prepared = req.prepare()
        # This header needs to be removed since most browsers do not allow us to set it manually
        del prepared.headers["Content-length"]
        # Perform the HTTP POST request
        s = requests.Session()
        return s.send(prepared)
