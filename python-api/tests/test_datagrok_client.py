from datetime import datetime
import unittest
from unittest.mock import patch, Mock, MagicMock
import pandas as pd
from io import StringIO
from datagrok_api import DatagrokClient, Group, GroupMembershipRequest

class TestDatagrokClient(unittest.TestCase):
    def setUp(self):
        self.client = DatagrokClient("Bearer token", "http://localhost:8082")
        self.mock_response = Mock()
        self.mock_response.json.return_value = {}
        self.mock_response.text = ""
        self.mock_response.content = b""

    @patch.object(DatagrokClient, '_request')
    def test_download_table(self, mock_request):
        # Mock CSV response
        csv_data = "col1,col2\n1,2\n3,4"
        self.mock_response.text = csv_data
        mock_request.return_value = self.mock_response

        df = self.client.download_table("test:table")
        
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 2)
        self.assertEqual(list(df.columns), ['col1', 'col2'])

    @patch.object(DatagrokClient, '_request')
    def test_upload_table(self, mock_request):
        # Mock response for successful upload
        self.mock_response.json.return_value = {"id": "table1", "name": "test:table"}
        mock_request.return_value = self.mock_response

        df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
        result = self.client.upload_table("test:table", df)
        
        self.assertEqual(result["id"], "table1")
        self.assertEqual(result["name"], "test:table")

    @patch.object(DatagrokClient, '_request')
    def test_download_file_csv(self, mock_request):
        # Mock CSV response
        csv_data = "col1,col2\n1,2\n3,4"
        self.mock_response.content = csv_data.encode('utf-8')
        mock_request.return_value = self.mock_response

        result = self.client.download_file("test:connector", "test.csv")
        
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)
        self.assertEqual(list(result.columns), ['col1', 'col2'])

    @patch.object(DatagrokClient, '_request')
    def test_download_file_binary(self, mock_request):
        # Mock binary response
        binary_data = b"binary content"
        self.mock_response.content = binary_data
        mock_request.return_value = self.mock_response

        result = self.client.download_file("test:connector", "test.bin")
        
        self.assertEqual(result, binary_data)

    @patch.object(DatagrokClient, '_request')
    def test_upload_file(self, mock_request):
        mock_request.return_value = self.mock_response

        with patch('builtins.open', unittest.mock.mock_open(read_data=b'test content')):
            self.client.upload_file("test:connector", "test.txt", "local.txt")
            
        mock_request.assert_called_once()

    @patch.object(DatagrokClient, '_request')
    def test_share_dashboard(self, mock_request):
        mock_request.return_value = self.mock_response

        self.client.share_dashboard("dashboard1", "group1,group2", "View")
        
        mock_request.assert_called_once()

    @patch.object(DatagrokClient, '_request')
    def test_create_dashboard(self, mock_request):
        self.mock_response.json.return_value = {"id": "dashboard1", "name": "test:dashboard"}
        mock_request.return_value = self.mock_response

        result = self.client.create_dashboard("test:dashboard", "table1,table2")
        
        self.assertEqual(result["id"], "dashboard1")
        self.assertEqual(result["name"], "test:dashboard")

    @patch.object(DatagrokClient, '_request')
    def test_call_function(self, mock_request):
        self.mock_response.json.return_value = {"result": "success"}
        mock_request.return_value = self.mock_response

        result = self.client.call_function("test:function", {"param1": "value1"})
        
        self.assertEqual(result, {"result": "success"})

    @patch.object(DatagrokClient, '_request')
    def test_list_groups(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = [
            {
                'id': 'group1',
                'name': 'Test Group 1',
                'friendlyName': 'Test Group 1',
                'description': 'Test Description 1',
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'personal': False,
                'hidden': False,
                'parents': [],
                'children': []
            },
            {
                'id': 'group2',
                'name': 'Test Group 2',
                'friendlyName': 'Test Group 2',
                'description': 'Test Description 2',
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'personal': False,
                'hidden': False,
                'parents': [],
                'children': []
            }
        ]
        mock_request.return_value = mock_response
        
        groups = self.client.list_groups(smart_filter="test", include_personal=True)
        
        self.assertEqual(len(groups), 2)
        self.assertEqual(groups[0].id, 'group1')
        self.assertEqual(groups[0].name, 'Test Group 1')
        self.assertEqual(groups[0].friendly_name, 'Test Group 1')
        self.assertEqual(groups[0].description, 'Test Description 1')
        self.assertEqual(groups[0].personal, False)
        self.assertEqual(groups[0].hidden, False)

    @patch.object(DatagrokClient, '_request')
    def test_lookup_groups(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = [
            {
                'id': 'group1',
                'name': 'Test Group',
                'friendlyName': 'Test Group',
                'description': 'Test Description',
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'personal': False,
                'hidden': False,
                'parents': [],
                'children': []
            }
        ]
        mock_request.return_value = mock_response
        
        groups = self.client.lookup_groups("Test Group")
        
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0].id, 'group1')
        self.assertEqual(groups[0].name, 'Test Group')
        self.assertEqual(groups[0].friendly_name, 'Test Group')
        self.assertEqual(groups[0].description, 'Test Description')

    @patch.object(DatagrokClient, '_request')
    def test_get_group(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = {
            'id': 'group1',
            'name': 'Test Group',
            'friendlyName': 'Test Group',
            'description': 'Test Description',
            'createdOn': '2025-04-29T09:53:40.502720Z',
            'updatedOn': '2025-04-29T09:53:40.502720Z',
            'personal': False,
            'hidden': False,
            'parents': [],
            'children': []
        }
        mock_request.return_value = mock_response
        
        group = self.client.get_group("group1")
        
        self.assertEqual(group.id, 'group1')
        self.assertEqual(group.name, 'Test Group')
        self.assertEqual(group.friendly_name, 'Test Group')
        self.assertEqual(group.description, 'Test Description')
        self.assertEqual(group.personal, False)
        self.assertEqual(group.hidden, False)

    @patch.object(DatagrokClient, '_request')
    def test_get_group_members(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = [
            {
                'id': 'user1',
                'name': 'User One',
                'friendlyName': 'User One',
                'description': 'Test User 1',
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'personal': True,
                'hidden': False,
                'parents': [],
                'children': []
            },
            {
                'id': 'user2',
                'name': 'User Two',
                'friendlyName': 'User Two',
                'description': 'Test User 2',
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'personal': True,
                'hidden': False,
                'parents': [],
                'children': []
            }
        ]
        mock_request.return_value = mock_response
        
        members = self.client.get_group_members("group1", admin=True)
        
        self.assertEqual(len(members), 2)
        self.assertEqual(members[0].id, 'user1')
        self.assertEqual(members[0].name, 'User One')
        self.assertEqual(members[0].personal, True)
        self.assertEqual(members[1].id, 'user2')
        self.assertEqual(members[1].name, 'User Two')
        self.assertEqual(members[1].personal, True)

    @patch.object(DatagrokClient, '_request')
    def test_get_group_memberships(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = [
            {
                'id': 'group1',
                'name': 'Test Group',
                'friendlyName': 'Test Group',
                'description': 'Test Description',
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'personal': False,
                'hidden': False,
                'parents': [],
                'children': []
            }
        ]
        mock_request.return_value = mock_response
        
        memberships = self.client.get_group_memberships("group1", admin=False)
        
        self.assertEqual(len(memberships), 1)
        self.assertEqual(memberships[0].id, 'group1')
        self.assertEqual(memberships[0].name, 'Test Group')
        self.assertEqual(memberships[0].personal, False)

    @patch.object(DatagrokClient, '_request')
    def test_request_group_membership(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = {
            'id': 'request1',
            'approved': None,
            'createdOn': '2025-04-29T09:53:40.502720Z',
            'updatedOn': '2025-04-29T09:53:40.502720Z',
            'resolutionDate': None,
            'from': {
                'id': 'user1',
                'name': 'User One',
                'friendlyName': 'User One',
                'description': 'Test User 1',
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'personal': True,
                'hidden': False,
                'parents': [],
                'children': []
            },
            'to': {
                'id': 'group1',
                'name': 'Test Group',
                'friendlyName': 'Test Group',
                'description': 'Test Description',
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'personal': False,
                'hidden': False,
                'parents': [],
                'children': []
            }
        }
        mock_request.return_value = mock_response
        
        request = self.client.request_group_membership("group1")
        
        self.assertIsInstance(request, GroupMembershipRequest)
        self.assertEqual(request.id, 'request1')
        self.assertIsNone(request.approved)
        self.assertIsNone(request.resolution_date)
        self.assertEqual(request.from_group.id, 'user1')
        self.assertEqual(request.to.id, 'group1')

    @patch.object(DatagrokClient, '_request')
    def test_get_group_membership_requests(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = [
            {
                'id': 'request1',
                'approved': None,
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'resolutionDate': None,
                'from': {
                    'id': 'user1',
                    'name': 'User One',
                    'friendlyName': 'User One',
                    'description': 'Test User 1',
                    'createdOn': '2025-04-29T09:53:40.502720Z',
                    'updatedOn': '2025-04-29T09:53:40.502720Z',
                    'personal': True,
                    'hidden': False,
                    'parents': [],
                    'children': []
                },
                'to': {
                    'id': 'group1',
                    'name': 'Test Group',
                    'friendlyName': 'Test Group',
                    'description': 'Test Description',
                    'createdOn': '2025-04-29T09:53:40.502720Z',
                    'updatedOn': '2025-04-29T09:53:40.502720Z',
                    'personal': False,
                    'hidden': False,
                    'parents': [],
                    'children': []
                }
            },
            {
                'id': 'request2',
                'approved': True,
                'createdOn': '2025-04-29T09:53:40.502720Z',
                'updatedOn': '2025-04-29T09:53:40.502720Z',
                'resolutionDate': '2025-04-29T10:00:00.000000Z',
                'from': {
                    'id': 'user2',
                    'name': 'User Two',
                    'friendlyName': 'User Two',
                    'description': 'Test User 2',
                    'createdOn': '2025-04-29T09:53:40.502720Z',
                    'updatedOn': '2025-04-29T09:53:40.502720Z',
                    'personal': True,
                    'hidden': False,
                    'parents': [],
                    'children': []
                },
                'to': {
                    'id': 'group1',
                    'name': 'Test Group',
                    'friendlyName': 'Test Group',
                    'description': 'Test Description',
                    'createdOn': '2025-04-29T09:53:40.502720Z',
                    'updatedOn': '2025-04-29T09:53:40.502720Z',
                    'personal': False,
                    'hidden': False,
                    'parents': [],
                    'children': []
                }
            }
        ]
        mock_request.return_value = mock_response
        
        requests = self.client.get_group_membership_requests("group1")
        
        self.assertEqual(len(requests), 2)
        self.assertEqual(requests[0].id, 'request1')
        self.assertIsNone(requests[0].approved)
        self.assertEqual(requests[1].id, 'request2')
        self.assertTrue(requests[1].approved)

    @patch.object(DatagrokClient, '_request')
    def test_approve_membership_request(self, mock_request):
        mock_request.return_value = self.mock_response
        
        self.client.approve_membership_request("request1")
        
        mock_request.assert_called_once()

    @patch.object(DatagrokClient, '_request')
    def test_deny_membership_request(self, mock_request):
        mock_request.return_value = self.mock_response
        
        self.client.deny_membership_request("request1")
        
        mock_request.assert_called_once()

    @patch.object(DatagrokClient, '_request')
    def test_get_current_user_group(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = {
            'id': 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6',
            'name': 'All users',
            'friendlyName': 'All users',
            'description': 'All users group',
            'createdOn': '2025-04-29T09:53:40.502720Z',
            'updatedOn': '2025-04-29T09:53:40.502720Z',
            'personal': False,
            'hidden': False,
            'parents': [],
            'children': []
        }
        mock_request.return_value = mock_response
        
        group = self.client.get_current_user_group()
        
        self.assertIsInstance(group, Group)
        self.assertEqual(group.id, 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6')
        self.assertEqual(group.name, 'All users')
        self.assertEqual(group.friendly_name, 'All users')
        self.assertEqual(group.description, 'All users group')
        self.assertEqual(group.personal, False)
        self.assertEqual(group.hidden, False)

    @patch.object(DatagrokClient, '_request')
    def test_save_group(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = {
            'id': 'group1',
            'name': 'Test Group',
            'friendlyName': 'Test Group',
            'description': 'Test Description',
            'createdOn': '2025-04-29T09:53:40.502720Z',
            'updatedOn': '2025-04-29T09:53:40.502720Z',
            'personal': False,
            'hidden': False,
            'parents': [],
            'children': []
        }
        mock_request.return_value = mock_response
        
        group = Group(self.client, id='group1', name='Test Group')
        saved_group = self.client.save_group(group, save_relations=True)
        
        self.assertEqual(saved_group.id, 'group1')
        self.assertEqual(saved_group.name, 'Test Group')
        self.assertEqual(saved_group.friendly_name, 'Test Group')
        self.assertEqual(saved_group.description, 'Test Description')
        self.assertEqual(saved_group.personal, False)
        self.assertEqual(saved_group.hidden, False)

if __name__ == '__main__':
    unittest.main()
